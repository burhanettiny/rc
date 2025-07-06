import streamlit as st
import re
import matplotlib.pyplot as plt
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction

st.set_page_config(page_title="PCR Primer Analizi", layout="wide")

st.title("🧬 DNA/RNA Primer ve Prob Analizi Aracı")

st.markdown("Bu araç, verdiğiniz DNA/RNA dizisi üzerinde primer ve prob yerleşimini gösterir, PCR bilgilerini hesaplar.")

# Kullanıcı girdileri
seq_input = st.text_area("🔢 DNA/RNA Sekansı", height=200).upper().replace(" ", "").replace("\n", "")
primer_f = st.text_input("➡️ Forward Primer").upper()
primer_r = st.text_input("⬅️ Reverse Primer").upper()
probe = st.text_input("🟨 Prob Dizisi (Opsiyonel)").upper()
molecule_type = st.selectbox("Molekül Tipi", ["DNA", "RNA"])

if molecule_type == "RNA":
    seq_input = seq_input.replace("U", "T")
    primer_f = primer_f.replace("U", "T")
    primer_r = primer_r.replace("U", "T")
    probe = probe.replace("U", "T")

def find_positions(seq, subseq):
    match = re.search(subseq, seq)
    return (match.start(), match.end()) if match else (-1, -1)

if st.button("🔍 Analizi Başlat"):

    if not seq_input or not primer_f or not primer_r:
        st.error("Sekans ve her iki primer zorunludur.")
    else:
        fwd_start, fwd_end = find_positions(seq_input, primer_f)
        rev_start, rev_end = find_positions(seq_input, primer_r)
        probe_start, probe_end = find_positions(seq_input, probe) if probe else (-1, -1)

        if fwd_start == -1 or rev_start == -1:
            st.error("Primer dizilerinden biri sekans içinde bulunamadı!")
        else:
            amplikon_uzunlugu = abs(rev_start - fwd_end)

            col1, col2 = st.columns(2)
            with col1:
                st.subheader("📐 Primer Bilgileri")
                st.write(f"**Forward Primer Pozisyonu:** {fwd_start}-{fwd_end}")
                st.write(f"**Reverse Primer Pozisyonu:** {rev_start}-{rev_end}")
                st.write(f"**Amplikon Uzunluğu:** {amplikon_uzunlugu} bp")
            with col2:
                st.subheader("🧪 PCR Koşulları")
                st.write(f"**Tm (Forward):** {mt.Tm_Wallace(primer_f):.2f} °C | GC: {gc_fraction(primer_f)*100:.2f}%")
                st.write(f"**Tm (Reverse):** {mt.Tm_Wallace(primer_r):.2f} °C | GC: {gc_fraction(primer_r)*100:.2f}%")
                if probe:
                    st.write(f"**Tm (Prob):** {mt.Tm_Wallace(probe):.2f} °C | GC: {gc(probe):.2f}%")

            # Görsel çizim
            fig, ax = plt.subplots(figsize=(12, 2))
            ax.plot(range(len(seq_input)), [1]*len(seq_input), 'k-', lw=1)

            ax.plot(range(fwd_start, fwd_end), [1.2]*len(primer_f), 'g>', lw=4, label='Forward Primer')
            ax.plot(range(rev_start, rev_end), [0.8]*len(primer_r), 'r<', lw=4, label='Reverse Primer')
            if probe and probe_start != -1:
                ax.plot(range(probe_start, probe_end), [1.4]*len(probe), 'b-', lw=4, label='Probe')

            ax.set_yticks([])
            ax.set_xlabel("Baz Pozisyonu")
            ax.set_title("Sekans Üzerinde Primer ve Prob Görseli")
            ax.legend()
            st.pyplot(fig)
