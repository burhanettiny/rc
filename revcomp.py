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
        highlight_regions = []

        num_blocks = len(seq_input) // 200 + 1
        for block in range(num_blocks):
            block_start = block * 200
            block_end = min((block + 1) * 200, len(seq_input))
            block_seq = seq_input[block_start:block_end]

            fig, ax = plt.subplots(figsize=(12, 2))
            ax.set_xlim(block_start, block_end)
            ax.set_ylim(0.6, 1.6)
            ax.set_yticks([])
            ax.set_xlabel("Baz Pozisyonu")
            ax.set_title(f"Primer ve Prob Yerleşimi (Pozisyon {block_start} - {block_end})")

            ax.plot(range(block_start, block_end), [1]*(block_end - block_start), 'k-', lw=1)

            for idx, s in enumerate(primer_sets):
                fwd = s["forward"]
                rev = s["reverse"]
                probe = s["probe"]
                rev_rc = reverse_complement(rev)

                fwd_start, fwd_end = find_positions(seq_input, fwd)
                rev_start, rev_end = find_positions(seq_input, rev_rc)
                probe_start, probe_end = find_positions(seq_input, probe) if probe else (-1, -1)

                if fwd_start != -1:
                    highlight_regions.append((fwd_start, fwd_end, primer_colors[idx]))
                if rev_start != -1:
                    highlight_regions.append((rev_start, rev_end, primer_colors[idx]))
                if probe and probe_start != -1:
                    highlight_regions.append((probe_start, probe_end, probe_colors[idx]))

                # Bilgileri ekle
                if fwd_start >= block_start and fwd_end <= block_end:
                    ax.plot(range(fwd_start, fwd_end), [1.2]*len(fwd), color=primer_colors[idx], lw=4, label=f"Set {idx+1} Forward")
                if rev_start >= block_start and rev_end <= block_end:
                    ax.plot(range(rev_start, rev_end), [0.8]*len(rev_rc), color=primer_colors[idx], lw=4, label=f"Set {idx+1} Reverse")
                if probe and probe_start >= block_start and probe_end <= block_end:
                    ax.plot(range(probe_start, probe_end), [1.4]*len(probe), color=probe_colors[idx], lw=4, label=f"Set {idx+1} Probe")

            # Harf bazlı sekans: her baz, kendi pozisyonuna göre renkli yazılır
            for i in range(block_start, block_end):
                base = seq_input[i]
                color = None
                for start, end, c in highlight_regions:
                    if start <= i < end:
                        color = c
                        break
                ax.text(i, 1.02, base, fontsize=7, color=color if color else 'black', ha='center', va='bottom', family='monospace')

            ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.25), ncol=3)
            st.pyplot(fig)
