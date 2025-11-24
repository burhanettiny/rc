import streamlit as st
from pypdf import PdfMerger, PdfReader, PdfWriter
from io import BytesIO
from docx import Document
import tempfile
import os

st.set_page_config(page_title="Belge Birleştirici", page_icon="📎", layout="centered")
st.title("📎 PDF & Word Birleştirici - Streamlit")
st.markdown("Bu uygulama PDF ve Word (DOCX) dosyalarını yükleyip sürükle-bırak yöntemiyle sırasını belirleyerek tek bir dosya haline getirir.")

st.markdown("---")

uploaded_files = st.file_uploader(
    "PDF veya Word dosyalarını yükleyin (çoklu seçim mümkün)",
    type=["pdf", "docx"],
    accept_multiple_files=True
)

# =====================================================================
#                       DOSYA SIRALAMA
# =====================================================================
if uploaded_files:
    file_names = [f.name for f in uploaded_files]
    st.subheader("Dosya sırası (sürükleyerek değiştirin)")

    sorted_files = st.sortable_items(file_names, key="file_sort")

    # =================================================================
    #                    PDF SAYFA YÖNETİMİ
    # =================================================================
    st.subheader("📄 PDF Sayfa Yönetimi")
    pdf_manage = st.selectbox("Sayfa yönetimi için bir PDF seçin", [n for n in file_names if n.lower().endswith('.pdf')])

    if pdf_manage:
        pdf_file = uploaded_files[file_names.index(pdf_manage)]
        reader = PdfReader(pdf_file)
        total_pages = len(reader.pages)

        st.write(f"Toplam sayfa: {total_pages}")

        page_list = [f"Sayfa {i+1}" for i in range(total_pages)]

        st.write("Sayfaları sürükleyerek yeniden sıralayın veya seçerek silin.")
        reordered = st.sortable_items(page_list, key=f"sort_pages_{pdf_manage}")
        delete_pages = st.multiselect("Silinecek sayfalar", reordered)

        if st.button("📌 Yeni PDF Üret (Sayfa Silme / Taşıma)"):
            writer = PdfWriter()
            for page_name in reordered:
                idx = int(page_name.split()[1]) - 1
                if page_name not in delete_pages:
                    writer.add_page(reader.pages[idx])

            out_pdf = BytesIO()
            writer.write(out_pdf)
            out_pdf.seek(0)

            st.success("Yeni PDF oluşturuldu!")
            st.download_button(
                "📥 Düzenlenmiş PDF'i İndir",
                out_pdf,
                f"edited_{pdf_manage}",
                mime="application/pdf",
            )

    # =================================================================
    #                        PDF BİRLEŞTİRME
    # =================================================================
    if st.button("🔀 PDF'leri Birleştir"):
        try:
            pdf_files = [
                uploaded_files[file_names.index(name)]
                for name in sorted_files if name.lower().endswith(".pdf")
            ]

            if not pdf_files:
                st.error("Birleştirilecek PDF dosyası bulunamadı.")
            else:
                merger = PdfMerger()
                for file in pdf_files:
                    merger.append(file)

                out = BytesIO()
                merger.write(out)
                merger.close()
                out.seek(0)

                st.success("PDF başarıyla birleştirildi!")
                st.download_button(
                    label="📥 Birleşmiş PDF'i İndir",
                    data=out,
                    file_name="merged.pdf",
                    mime="application/pdf",
                )
        except Exception as e:
            st.error(f"PDF birleştirme hatası: {e}")

    # =================================================================
    #                      WORD BİRLEŞTİRME (DOCX)
    # =================================================================
    if st.button("📝 Word (DOCX) Birleştir"):
        try:
            word_files = [
                uploaded_files[file_names.index(name)]
                for name in sorted_files if name.lower().endswith(".docx")
            ]

            if not word_files:
                st.error("Birleştirilecek Word dosyası bulunamadı.")
            else:
                merged_doc = Document()
                first = True

                for file in word_files:
                    temp_path = tempfile.mktemp(suffix=".docx")
                    with open(temp_path, "wb") as tmp:
                        tmp.write(file.getbuffer())
                    sub_doc = Document(temp_path)

                    if first:
                        for p in sub_doc.paragraphs:
                            merged_doc.add_paragraph(p.text)
                        first = False
                    else:
                        merged_doc.add_page_break()
                        for p in sub_doc.paragraphs:
                            merged_doc.add_paragraph(p.text)

                    os.remove(temp_path)

                out_path = tempfile.mktemp(suffix=".docx")
                merged_doc.save(out_path)

                with open(out_path, "rb") as f:
                    st.success("Word belgeleri birleştirildi!")
                    st.download_button(
                        "📥 Birleşmiş Word Belgesini İndir",
                        f.read(),
                        "merged.docx",
                        mime="application/vnd.openxmlformats-officedocument.wordprocessingml.document",
                    )
        except Exception as e:
            st.error(f"Word birleştirme hatası: {e}")

    # =================================================================
    #                DOCX + PDF → TEK PDF BİRLEŞTİRME
    # =================================================================
    if st.button("📄 DOCX + PDF → Tek PDF Birleştir"):
        try:
            import docx2pdf

            pdf_files = [name for name in sorted_files if name.lower().endswith(".pdf")]
            docx_files = [name for name in sorted_files if name.lower().endswith(".docx")]

            temp_pdf_list = []

            # DOCX → PDF dönüşümü
            for name in docx_files:
                file = uploaded_files[file_names.index(name)]
                tmp_docx = tempfile.mktemp(suffix=".docx")
                with open(tmp_docx, "wb") as tmp:
                    tmp.write(file.getbuffer())
                tmp_pdf = tempfile.mktemp(suffix=".pdf")
                docx2pdf.convert(tmp_docx, tmp_pdf)
                temp_pdf_list.append(tmp_pdf)

            # PDF birleştirme
            merger = PdfMerger()

            # Sıraya göre ekle
            temp_idx = 0
            for name in sorted_files:
                if name.lower().endswith(".pdf"):
                    merger.append(uploaded_files[file_names.index(name)])
                else:
                    merger.append(temp_pdf_list[temp_idx])
                    temp_idx += 1

            out = BytesIO()
            merger.write(out)
            out.seek(0)
            merger.close()

            st.success("DOCX + PDF birlikte tek PDF olarak birleştirildi!")
            st.download_button(
                "📥 Tek PDF Olarak İndir",
                out,
                "merged_all.pdf",
                mime="application/pdf",
            )
        except Exception as e:
            st.error(f"Birleştirme hatası: {e}")

else:
    st.info("Başlamak için PDF veya Word dosyalarını yükleyin.")

st.markdown("---")
st.caption("Not: Çok büyük dosyalarda bellek sınırları sorun oluşturabilir. Yerel çalıştırma daha stabil olabilir.")

st.markdown("""
**Gereksinimler**:
- `pip install streamlit`
- `pip install pypdf`
- `pip install python-docx`
- `pip install docx2pdf`
- `pip install streamlit-sortable`

**Çalıştırma**:
