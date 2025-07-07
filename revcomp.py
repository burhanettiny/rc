import streamlit as st
import re
import matplotlib.pyplot as plt
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction

st.set_page_config(page_title="Multiplex PCR Primer Analizi", layout="wide")
st.title("🧬 Multiplex DNA/RNA Primer ve Prob Analizi Aracı")

st.markdown("Bu araç, birden fazla primer/prob seti ile PCR analizi yapmanızı sağlar. Tm, GC, Ta, amplikon uzunluğu hesaplanır ve görsel olarak gösterilir.")

# Kullanıcı girdileri
seq_input = st.text_area("🔢 DNA/RNA Sekansı", height=200).upper().replace(" ", "").replace("\n", "")
molecule_type = st.selectbox("Molekül Tipi", ["DNA", "RNA"])

if molecule_type == "RNA":
    seq_input = seq_input.replace("U", "T")

# Multiplex giriş
primer_set_count = st.number_input("🔢 Primer Seti Sayısı (Multiplex)", min_value=1, max_value=5, value=1, step=1)
primer_sets = []

for i in range(primer_set_count):
    with st.expander(f"🧬 Primer Set {i+1}"):
        fwd = st.text_input(f"➡️ Forward Primer {i+1} (5'-3')", key=f"fwd_{i}").upper().replace("U", "T")
        rev = st.text_input(f"⬅️ Reverse Primer {i+1} (5'-3')", key=f"rev_{i}").upper().replace("U", "T")
        prb = st.text_input(f"🟨 Prob {i+1} (opsiyonel)", key=f"prb_{i}").upper().replace("U", "T")
        primer_sets.append({"forward": fwd, "reverse": rev, "probe": prb})

def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

def find_positions(seq, subseq):
    match = re.search(subseq, seq)
    return (match.start(), match.end()) if match else (-1, -1)

primer_colors = ['green', 'red', 'purple', 'orange', 'brown']
probe_colors = ['cyan', 'magenta', 'gray', 'blue', 'olive']

if st.button("🔍 Analizi Başlat"):
    if not seq_input:
        st.error("Sekans girişi zorunludur.")
    elif any(not s["forward"] or not s["reverse"] for s in primer_sets):
        st.error("Her primer seti için forward ve reverse primer gereklidir.")
    else:
        st.subheader("📐 Primer Seti Analizi")
        fig, ax = plt.subplots(figsize=(12, 2))
        ax.plot(range(len(seq_input)), [1]*len(seq_input), 'k-', lw=1)

        for idx, s in enumerate(primer_sets):
            fwd = s["forward"]
            rev = s["reverse"]
            probe = s["probe"]

            rev_rc = reverse_complement(rev)
            st.info(f"Set {idx+1}: Reverse primer 5'-3' yönünde girildiği için reverse complement alındı.")

            fwd_start, fwd_end = find_positions(seq_input, fwd)
            rev_start, rev_end = find_positions(seq_input, rev_rc)
            probe_start, probe_end = find_positions(seq_input, probe) if probe else (-1, -1)

            if fwd_start == -1 or rev_start == -1:
                st.error(f"Set {idx+1}: Primerlar sekans içinde bulunamadı!")
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

            ax.plot(range(fwd_start, fwd_end), [1.2]*len(fwd), color=primer_colors[idx], lw=4, label=f"Set {idx+1} Forward")
            ax.plot(range(rev_start, rev_end), [0.8]*len(rev_rc), color=primer_colors[idx], lw=4, label=f"Set {idx+1} Reverse")
            if probe and probe_start != -1:
                ax.plot(range(probe_start, probe_end), [1.4]*len(probe), color=probe_colors[idx], lw=4, label=f"Set {idx+1} Probe")

        ax.set_yticks([])
        ax.set_xlabel("Baz Pozisyonu")
        ax.set_title("Sekans Üzerinde Primer ve Prob Yerleşimi")
        ax.legend()
        st.pyplot(fig)

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

        st.pyplot
