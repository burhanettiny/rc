import streamlit as st
from pypdf import PdfMerger
from io import BytesIO

st.set_page_config(page_title="PDF Birleştirici", page_icon="📎", layout="centered")
st.title("📎 PDF Birleştirici - Streamlit")
st.markdown("Bu uygulama birden fazla PDF dosyasını yükleyip sırasını belirleyerek tek bir PDF haline getirir.")

st.markdown("---")

uploaded_files = st.file_uploader("PDF dosyalarını yükleyin (çoklu seçim mümkün)", type=['pdf'], accept_multiple_files=True)

if uploaded_files:
    # Dosyaları listele
    st.subheader("Yüklenen dosyalar")
    file_names = [f.name for f in uploaded_files]
    for i, name in enumerate(file_names, start=1):
        st.write(f"{i}. {name}")

    st.info("Sırayı değiştirmek için dosya indekslerini virgülle ayırarak (ör. 2,1,3) girin. Varsayılan sıra yükleme sırasıdır.")
    order_input = st.text_input("Birleştirme sırası (indekslerle)", value=",")

    # Varsayılan sıra: 1,2,3... gösterimi
    default_order = ",".join(str(i) for i in range(1, len(uploaded_files) + 1))
    if order_input.strip() == ",":
        order_input = default_order

    try:
        # Parse order
        indices = [int(x.strip()) for x in order_input.split(",") if x.strip()]
        if sorted(indices) != list(range(1, len(uploaded_files) + 1)):
            st.warning("Girdiğiniz sıra bütün dosya indekslerini içermiyor veya tekrar içeriyor. Lütfen geçerli bir sıra girin.")
        else:
            # Butona basılınca birleştir
            if st.button("🔀 PDFleri Birleştir"):
                merger = PdfMerger()
                try:
                    for idx in indices:
                        file_obj = uploaded_files[idx - 1]
                        # file_obj: UploadedFile
                        merger.append(file_obj)

                    out = BytesIO()
                    merger.write(out)
                    merger.close()
                    out.seek(0)

                    st.success("PDF başarıyla birleştirildi!")
                    merged_name = "merged_" + "_".join([uploaded_files[i-1].name.replace(' ', '_') for i in indices])
                    if not merged_name.lower().endswith('.pdf'):
                        merged_name += '.pdf'

                    st.download_button(label="📥 Birleşmiş PDF'i İndir", data=out, file_name=merged_name, mime='application/pdf')
                except Exception as e:
                    st.error(f"Birleştirme sırasında hata oluştu: {e}")
    except ValueError:
        st.error("Sıra girdisini okurken hata: Lütfen sadece virgülle ayrılmış sayılar girin (ör. 1,2,3).")

else:
    st.info("Başlamak için soldan veya yukarıdan PDF dosyaları yükleyin.")

st.markdown("---")
st.caption("Not: Sunucuda çok büyük dosyalar yüklenmesi bellek/kota sorunlarına yol açabilir. Yerel çalıştırmada daha yüksek limitler için Streamlit ayarlarınıza bakın.")

# Yardım / Gereksinimler
st.markdown("**Gereksinimler**: `pip install streamlit pypdf`\n\n**Çalıştırma**: `streamlit run streamlit_pdf_birlestir.py`")
