import streamlit as st
import re
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction

# 📌 IUPAC kodları
IUPAC_CODES = {
    'A': 'A', 'T': 'T', 'G': 'G', 'C': 'C',
    'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
    'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
    'H': '[ACT]', 'V': '[ACG]', 'N': '[ATGC]'
}

# 🔎 Motifi regex'e çevir
def iupac_to_regex(motif):
    return ''.join(IUPAC_CODES.get(base, base) for base in motif)

# 🔬 Enzim kesim bölgesi bul
def find_all_enzymes_in_sequence(seq, enzymes_dict):
    results = {}
    for enzyme, motif in enzymes_dict.items():
        pattern = iupac_to_regex(motif)
        matches = list(re.finditer(pattern, seq))
        if matches:
            results[enzyme] = [(m.start(), m.end()) for m in matches]
    return results

# 📌 Restriksiyon enzimleri
RE_SITES = {
    "EcoRI": "GAATTC", "BamHI": "GGATCC", "HindIII": "AAGCTT", "NotI": "GCGGCCGC",
    "XhoI": "CTCGAG", "PstI": "CTGCAG", "SacI": "GAGCTC", "SalI": "GTCGAC",
    "SmaI": "CCCGGG", "KpnI": "GGTACC", "ApaI": "GGGCCC", "NcoI": "CCATGG",
    "MluI": "ACGCGT", "NheI": "GCTAGC", "SpeI": "ACTAGT", "ClaI": "ATCGAT",
    "DraI": "TTTAAA", "BglII": "AGATCT", "XbaI": "TCTAGA", "EagI": "CGGCCG",
    "TestEnz": "CCANNNNNNTGG"
}

def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

def find_positions(seq, subseq):
    match = re.search(subseq, seq)
    return (match.start(), match.end()) if match else (-1, -1)

def find_methylation_regions(seq, motif, gap_threshold=5):
    matches = [m.start() for m in re.finditer(f'(?={motif})', seq)]
    if not matches:
        return []
    regions = []
    group = [matches[0]]
    for pos in matches[1:]:
        if pos - group[-1] <= gap_threshold:
            group.append(pos)
        else:
            start, end = group[0], group[-1] + len(motif)
            count = len(group)
            percent = min(100, count * 20)
            regions.append({"start": start, "end": end, "count": count, "percent": percent})
            group = [pos]
    start, end = group[0], group[-1] + len(motif)
    count = len(group)
    percent = min(100, count * 20)
    regions.append({"start": start, "end": end, "count": count, "percent": percent})
    return regions

def highlight_sequence(seq, primer_sets, methylation_regions=None, enzyme_sites=None, line_length=80):
    style = """
    <style>
    .seq-box { font-family: Courier New, monospace; font-size: 14px; line-height: 1.4; white-space: pre-wrap; }
    .fwd { background-color: #90ee90; }
    .rev { background-color: #ffcccb; }
    .prb { background-color: #add8e6; }
    .met { background-color: #ffff99; }
    .enz { background-color: #dda0dd; }
    </style>
    """
    seq_list = list(seq)
    tags = [''] * len(seq)
    complement_map = str.maketrans("ATGC", "TACG")
    comp_seq = seq.translate(complement_map)
    comp_list = list(comp_seq)
    for s in primer_sets:
        for label, primer in [('fwd', s["forward"]), ('rev', s["reverse"]), ('prb', s["probe"])]:
            if not primer:
                continue
            rev_comp = reverse_complement(primer)
            if primer in seq:
                start = seq.find(primer)
                for i in range(start, start + len(primer)):
                    tags[i] = label
            elif rev_comp in seq:
                start = seq.find(rev_comp)
                for i in range(start, start + len(primer)):
                    tags[i] = label
    if methylation_regions:
        for region in methylation_regions:
            for i in range(region["start"], region["end"]):
                if tags[i] == '':
                    tags[i] = 'met'
    if enzyme_sites:
        for positions in enzyme_sites.values():
            for start, end in positions:
                for i in range(start, end):
                    if tags[i] == '':
                        tags[i] = 'enz'
    def render_line(bases, tags):
        return ''.join(f'<span class="{tag}">{base}</span>' if tag else base for base, tag in zip(bases, tags))
    lines = []
    for i in range(0, len(seq), line_length):
        top_line = render_line(seq_list[i:i+line_length], tags[i:i+line_length])
        bottom_line = render_line(comp_list[i:i+line_length], ['']*min(line_length,len(seq)-i))
        lines.append(f"5' {top_line} 3'<br>3' {bottom_line} 5'")
    return f"{style}<div class='seq-box'>" + "<br><br>".join(lines) + "</div>"

# 🧬 Streamlit Arayüz
st.set_page_config(page_title="Multiplex PCR Primer Analizi", layout="wide")
st.title("🧬 Sekansta Restriksiyon Enzim & Metilasyon & Primer Analizi")

seq_input = st.text_area("🔢 DNA/RNA Sekansı", height=200).upper().replace(" ", "").replace("\n", "")
molecule_type = st.selectbox("Molekül Tipi", ["DNA", "RNA"])
if molecule_type == "RNA":
    seq_input = seq_input.replace("U", "T")

methylation_motif = st.text_input("🧬 Metilasyon Motifi (örn: CG)", value="CG").upper().replace("U", "T")

primer_set_count = st.number_input("🔢 Primer Seti Sayısı", min_value=1, max_value=5, value=1, step=1)
primer_sets = []
for i in range(primer_set_count):
    with st.expander(f"🧬 Primer Set {i+1}"):
        fwd = st.text_input(f"➡️ Forward Primer {i+1}", key=f"fwd_{i}").upper().replace("U", "T")
        rev = st.text_input(f"⬅️ Reverse Primer {i+1}", key=f"rev_{i}").upper().replace("U", "T")
        prb = st.text_input(f"🟨 Prob {i+1} (opsiyonel)", key=f"prb_{i}").upper().replace("U", "T")
        primer_sets.append({"forward": fwd, "reverse": rev, "probe": prb})

if st.button("🔍 Analizi Başlat"):
    if not seq_input:
        st.error("⚠️ Sekans girişi zorunludur.")
    elif any(not s["forward"] or not s["reverse"] for s in primer_sets):
        st.error("⚠️ Her primer seti için forward ve reverse primer girilmelidir.")
    else:
        st.subheader("📐 Primer Analizi")
        for idx, s in enumerate(primer_sets):
            fwd = s["forward"]
            rev = s["reverse"]
            rev_rc = reverse_complement(rev)
            fwd_start, fwd_end = find_positions(seq_input, fwd)
            rev_start, rev_end = find_positions(seq_input, rev_rc)
            if fwd_start == -1 or rev_start == -1:
                st.warning(f"Set {idx+1}: Primerlar bulunamadı.")
                continue
            amplikon = rev_end - fwd_start
            Tm_f = mt.Tm_Wallace(fwd)
            Tm_r = mt.Tm_Wallace(rev)
            Ta = ((Tm_f + Tm_r) / 2) - 5
            st.markdown(f"### 🧪 Primer Set {idx+1}")
            st.write(f"**Forward Primer Pozisyonu:** {fwd_start}-{fwd_end}")
            st.write(f"**Reverse Primer Pozisyonu:** {rev_start}-{rev_end}")
            st.write(f"**Amplikon Uzunluğu:** {amplikon} bp")
            st.write(f"**Tm (Forward):** {Tm_f:.2f} °C | GC: {gc_fraction(fwd)*100:.2f}%")
            st.write(f"**Tm (Reverse):** {Tm_r:.2f} °C | GC: {gc_fraction(rev)*100:.2f}%")
            if s["probe"]:
                st.write(f"**Tm (Prob):** {mt.Tm_Wallace(s['probe']):.2f} °C | GC: {gc_fraction(s['probe'])*100:.2f}%")
            st.markdown(f"**🔹 Optimum Annealing Temperature (Ta):** {Ta:.2f} °C 🔥")

        enzyme_sites = find_all_enzymes_in_sequence(seq_input, RE_SITES)
        methylation_regions = find_methylation_regions(seq_input, methylation_motif)

        st.subheader("🧬 Renklendirilmiş Sekans Görünümü")
        html = highlight_sequence(seq_input, primer_sets, methylation_regions, enzyme_sites)
        st.markdown(html, unsafe_allow_html=True)

        st.subheader("🔪 Restriksiyon Enzim Kesim Bölgeleri")
        if enzyme_sites:
            for enz, positions in enzyme_sites.items():
                st.write(f"🧬 {enz}: " + ", ".join([f"{start}-{end}" for start, end in positions]))
        else:
            st.info("Hiçbir enzim kesim bölgesi bulunamadı.")

        st.subheader("🧬 Metilasyon Bölgeleri")
        if methylation_regions:
            df = pd.DataFrame([{"Başlangıç": r["start"], "Bitiş": r["end"], "Motif Sayısı": r["count"], "% Metilasyon": r["percent"]} for r in methylation_regions])
            st.dataframe(df)
        else:
            st.info("Hiçbir metilasyon bölgesi bulunamadı.")
