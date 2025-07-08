import streamlit as st
import re
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction

# Restriksiyon enzim motifleri
RE_SITES = {
    "EcoRI": "GAATTC", "BamHI": "GGATCC", "HindIII": "AAGCTT", "NotI": "GCGGCCGC",
    "XhoI": "CTCGAG", "PstI": "CTGCAG", "SacI": "GAGCTC", "SalI": "GTCGAC",
    "SmaI": "CCCGGG", "KpnI": "GGTACC", "ApaI": "GGGCCC", "NcoI": "CCATGG",
    "MluI": "ACGCGT", "NheI": "GCTAGC", "SpeI": "ACTAGT", "ClaI": "ATCGAT",
    "DraI": "TTTAAA", "BglII": "AGATCT", "XbaI": "TCTAGA", "EagI": "CGGCCG",
    "SbfI": "CCTGCAGG", "AscI": "GGCGCGCC", "PacI": "TTAATTAA", "SwaI": "ATTTAAAT",
    "MfeI": "CAATTG", "NsiI": "ATGCAT", "AflII": "CTTAAG", "BspEI": "TCCGGA",
    "StuI": "AGGCCT", "ScaI": "AGTACT"
}

st.set_page_config(page_title="Multiplex PCR Primer Analizi", layout="wide")
st.title("🧬 Sekansta Otomatik Restriksiyon Enzim ve Metilasyon Analizi")

seq_input = st.text_area("🔢 DNA/RNA Sekansı", height=200).replace(" ", "").replace("\n", "")
molecule_type = st.selectbox("Molekül Tipi", ["DNA", "RNA"])
if molecule_type == "RNA":
    seq_input = seq_input.replace("U", "T")

seq_input_upper = seq_input.upper()

methylation_motif = st.text_input("🧬 Metilasyon Motifi veya Sekansı (örn: CG)", value="CG").upper().replace("U", "T")

primer_set_count = st.number_input("🔢 Primer Seti Sayısı (Multiplex)", min_value=1, max_value=5, value=1, step=1)
primer_sets = []
for i in range(primer_set_count):
    with st.expander(f"🧬 Primer Set {i+1}"):
        fwd = st.text_input(f"➡️ Forward Primer {i+1} (5'-3')", key=f"fwd_{i}").replace("U", "T").upper()
        rev = st.text_input(f"⬅️ Reverse Primer {i+1} (5'-3')", key=f"rev_{i}").replace("U", "T").upper()
        prb = st.text_input(f"🟨 Prob {i+1} (opsiyonel)", key=f"prb_{i}").replace("U", "T").upper()
        primer_sets.append({"forward": fwd, "reverse": rev, "probe": prb})

def reverse_complement(seq):
    return str(Seq(seq).reverse_complement()).upper()

def find_positions(seq, subseq):
    match = re.search(subseq, seq, re.IGNORECASE)
    return (match.start(), match.end()) if match else (-1, -1)

def find_all_enzymes_in_sequence(seq, enzymes_dict):
    results = {}
    for enzyme, motif in enzymes_dict.items():
        matches = list(re.finditer(motif, seq, re.IGNORECASE))
        if matches:
            results[enzyme] = [(m.start(), m.end()) for m in matches]
    return results

def find_methylation_regions(seq, motif, gap_threshold=5):
    matches = [m.start() for m in re.finditer(f'(?={motif})', seq, re.IGNORECASE)]
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
    seq_upper = seq.upper()
    seq_list = list(seq_upper)
    tags = [''] * len(seq_upper)
    complement_map = str.maketrans("ATGC", "TACG")
    comp_seq = seq_upper.translate(complement_map)
    comp_list = list(comp_seq)

    for s in primer_sets:
        for label, primer in [('fwd', s["forward"]), ('rev', s["reverse"]), ('prb', s["probe"])]:
            if not primer:
                continue
            rev_comp = reverse_complement(primer)
            primer_lower = primer.lower()
            rev_comp_lower = rev_comp.lower()
            seq_lower = seq.lower()
            if primer_lower in seq_lower:
                start = seq_lower.find(primer_lower)
                for i in range(start, start + len(primer)):
                    tags[i] = label
            elif rev_comp_lower in seq_lower:
                start = seq_lower.find(rev_comp_lower)
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
    for i in range(0, len(seq_upper), line_length):
        top_line = render_line(seq_list[i:i+line_length], tags[i:i+line_length])
        bottom_line = render_line(comp_list[i:i+line_length], ['']*min(line_length,len(seq_upper)-i))
        lines.append(f"5' {top_line} 3'<br>3' {bottom_line} 5'")

    return f"{style}<div class='seq-box'>" + "<br><br>".join(lines) + "</div>"

if st.button("🔍 Analizi Başlat"):
    if not seq_input:
        st.error("Sekans girişi zorunludur.")
    elif any(not s["forward"] or not s["reverse"] for s in primer_sets):
        st.error("Her primer seti için forward ve reverse primer gereklidir.")
    else:
        st.subheader("📐 Primer Seti Analizi")
        for idx, s in enumerate(primer_sets):
            fwd = s["forward"]
            rev = s["reverse"]
            probe = s["probe"]
            rev_rc = reverse_complement(rev)
            fwd_start, fwd_end = find_positions(seq_input_upper, fwd)
            rev_start, rev_end = find_positions(seq_input_upper, rev_rc)
            if fwd_start == -1 or rev_start == -1:
                st.warning(f"Set {idx+1}: Primerlar sekans içinde bulunamadı.")
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
            if probe:
                st.write(f"**Tm (Prob):** {mt.Tm_Wallace(probe):.2f} °C | GC: {gc_fraction(probe)*100:.2f}%")
            st.markdown(f"**🔹 Optimum Annealing Temperature (Ta):** {Ta:.2f} °C 🔥")

        enzyme_sites = find_all_enzymes_in_sequence(seq_input_upper, RE_SITES)
        methylation_regions = find_methylation_regions(seq_input_upper, methylation_motif)

        st.subheader("🧬 Sekans Görünümü")
        html = highlight_sequence(seq_input_upper, primer_sets, methylation_regions, enzyme_sites)
        st.markdown(html, unsafe_allow_html=True)

        st.subheader("🔪 Restriksiyon Enzim Kesim Bölgeleri")
        if enzyme_sites:
            for enzyme, positions in enzyme_sites.items():
                st.write(f"🧬 {enzyme}:")
                for i, (start, end) in enumerate(positions, 1):
                    st.markdown(f"- Kesim bölgesi {i}: {start} - {end}")
        else:
            st.info("Sekans içinde bilinen enzim kesim bölgesi bulunamadı.")

        st.subheader("🧬 Metilasyon Bölgeleri (Otomatik)")
        if methylation_regions:
            df = pd.DataFrame([{"Başlangıç": r["start"], "Bitiş": r["end"], "Motif Sayısı": r["count"], "% Metilasyon": r["percent"]} for r in methylation_regions])
            st.dataframe(df)
        else:
            st.info("Belirtilen motife göre metilasyon bölgesi bulunamadı.")
