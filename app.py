import streamlit as st
import re
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction

# Sayfa ayarları
st.set_page_config(page_title="Multiplex PCR Primer Analizi", layout="wide")
st.title("🧬 Multiplex DNA/RNA Primer ve Prob Analizi Aracı")

st.markdown("Bu araç, birden fazla primer/prob seti ile PCR analizi yapmanızı sağlar. Tm, GC, Ta, amplikon uzunluğu hesaplanır ve sekans üzerinde renkli gösterilir.")

# Kullanıcı girişi
seq_input = st.text_area("🔢 DNA/RNA Sekansı", height=200).upper().replace(" ", "").replace("\n", "")
molecule_type = st.selectbox("Molekül Tipi", ["DNA", "RNA"])

if molecule_type == "RNA":
    seq_input = seq_input.replace("U", "T")

# Primer set sayısı
primer_set_count = st.number_input("🔢 Primer Seti Sayısı (Multiplex)", min_value=1, max_value=5, value=1, step=1)
primer_sets = []

# Kullanıcıdan primer/prob girişleri
for i in range(primer_set_count):
    with st.expander(f"🧬 Primer Set {i+1}"):
        fwd = st.text_input(f"➡️ Forward Primer {i+1} (5'-3')", key=f"fwd_{i}").upper().replace("U", "T")
        rev = st.text_input(f"⬅️ Reverse Primer {i+1} (5'-3')", key=f"rev_{i}").upper().replace("U", "T")
        prb = st.text_input(f"🟨 Prob {i+1} (opsiyonel)", key=f"prb_{i}").upper().replace("U", "T")
        primer_sets.append({"forward": fwd, "reverse": rev, "probe": prb})

# Yardımcı fonksiyonlar
def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

def find_positions(seq, subseq):
    match = re.search(subseq, seq)
    return (match.start(), match.end()) if match else (-1, -1)

# Sekans vurgulayıcı HTML görsel fonksiyonu
def highlight_sequence(seq, primer_sets, molecule_type, line_length=200):
    style = """
    <style>
    .seq-box { font-family: Courier New, monospace; font-size: 15px; line-height: 1.6; white-space: pre-wrap; }
    .fwd { background-color: #90ee90; }
    .rev { background-color: #ffcccb; }
    .prb { background-color: #add8e6; }
    </style>
    """

    seq_list = list(seq)
    top_tags = [''] * len(seq)
    bottom_tags = [''] * len(seq)

    # Komplement dizi (alt zincir)
    complement_map = str.maketrans("ATGC", "TACG")
    comp_seq = seq.translate(complement_map)
    comp_list = list(comp_seq)

    for s in primer_sets:
        for label, primer in [('fwd', s["forward"]), ('rev', s["reverse"]), ('prb', s["probe"])]:
            if not primer:
                continue

            rev_comp = reverse_complement(primer)
            # Primer doğrudan referans dizideyse → üst zincire boya
            if primer in seq:
                start = seq.find(primer)
                for i in range(start, start + len(primer)):
                    top_tags[i] = label
            # Reverse complement olarak varsa → alt zincire boya
            elif rev_comp in seq:
                start = seq.find(rev_comp)
                for i in range(start, start + len(primer)):
                    bottom_tags[i] = label
            else:
                continue  # eşleşme yoksa boyama yok

    def render_line(bases, tags):
        return ''.join(
            f'<span class="{tag}">{base}</span>' if tag else base
            for base, tag in zip(bases, tags)
        )

    # Satır satır gösterim
    lines = []
    for i in range(0, len(seq), line_length):
        line_seq = seq_list[i:i+line_length]
        line_top_tags = top_tags[i:i+line_length]

        line_comp = comp_list[i:i+line_length]
        line_bottom_tags = bottom_tags[i:i+line_length]

        top_line = render_line(line_seq, line_top_tags)
        bottom_line = render_line(line_comp, line_bottom_tags)

        lines.append(f"5' {top_line} 3'<br>3' {bottom_line} 5'")

    html = f"{style}<div class='seq-box'>" + "<br><br>".join(lines) + "</div>"
    return html

# Analiz başlatma
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

            fwd_start, fwd_end = find_positions(seq_input, fwd)
            rev_start, rev_end = find_positions(seq_input, rev_rc)
            probe_start, probe_end = find_positions(seq_input, probe) if probe else (-1, -1)

            if fwd_start == -1 or rev_start == -1:
                st.warning(f"Set {idx+1}: Primerlar sekans içinde bulunamadı.")
                continue

            amplikon = rev_end - fwd_start
            Tm_f = mt.Tm_Wallace(fwd)
            Tm_r = mt.Tm_Wallace(rev)
            Ta = ((Tm_f + Tm_r) / 2) - 5

            with st.container():
                st.markdown(f"### 🧪 Primer Set {idx+1}")
                st.write(f"**Forward Primer Pozisyonu:** {fwd_start}-{fwd_end}")
                st.write(f"**Reverse Primer Pozisyonu:** {rev_start}-{rev_end}")
                st.write(f"**Amplikon Uzunluğu:** {amplikon} bp")
                st.write(f"**Tm (Forward):** {Tm_f:.2f} °C | GC: {gc_fraction(fwd)*100:.2f}%")
                st.write(f"**Tm (Reverse):** {Tm_r:.2f} °C | GC: {gc_fraction(rev)*100:.2f}%")
                if probe:
                    st.write(f"**Tm (Prob):** {mt.Tm_Wallace(probe):.2f} °C | GC: {gc_fraction(probe)*100:.2f}%")
                st.markdown(f"**🔹 Optimum Annealing Temperature (Ta):** {Ta:.2f} °C 🔥")
                st.markdown("> **Formül:**  \n> Ta = ((Tm_forward + Tm_reverse) / 2) - 5")

        # Sekans görseli (yeni sistem)
        st.subheader("🧬 Sekans Üzerinde Primer ve Prob Yerleşimi")
        html = highlight_sequence(seq_input, primer_sets, molecule_type)
        st.markdown(html, unsafe_allow_html=True)

        # PCR Döngüsü
        st.subheader("🎛️ PCR Döngüsü Özelleştir")
        col1, col2 = st.columns(2)
        with col1:
            denaturation_temp = st.number_input("Denatürasyon Sıcaklığı (°C)", value=95)
            denaturation_time = st.number_input("Denatürasyon Süresi (sn)", value=30)
        with col2:
            annealing_time = st.number_input("Annealing Süresi (sn)", value=30)
            extension_temp = st.number_input("Uzama Sıcaklığı (°C)", value=72)
            extension_time = st.number_input("Uzama Süresi (sn)", value=60)
        cycle_count = st.slider("🔁 Döngü Sayısı", min_value=10, max_value=50, value=35)

        Ta = ((Tm_f + Tm_r) / 2) - 5 if 'Tm_f' in locals() else 60
        pcr_table = pd.DataFrame({
            "Adım": ["Denatürasyon", "Annealing", "Uzama"],
            "Sıcaklık (°C)": [denaturation_temp, Ta, extension_temp],
            "Süre (sn)": [denaturation_time, annealing_time, extension_time]
        })

        st.subheader("📋 Özelleştirilmiş PCR Döngüsü")
        st.table(pcr_table)
