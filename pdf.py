import streamlit as st
import pandas as pd
import tabula
import tempfile
import os

st.set_page_config(page_title="PDF to Excel", layout="centered")
st.title("📄 PDF Tabloyu Excel'e Dönüştür")

uploaded_file = st.file_uploader("PDF dosyasını yükle", type=["pdf"])

if uploaded_file:
    with tempfile.NamedTemporaryFile(delete=False, suffix=".pdf") as tmp_file:
        tmp_file.write(uploaded_file.read())
        tmp_pdf_path = tmp_file.name

    try:
        st.info("PDF içindeki tablolar okunuyor...")
        tables = tabula.read_pdf(tmp_pdf_path, pages='all', multiple_tables=True)

        if not tables:
            st.warning("PDF içinde tablo bulunamadı.")
        else:
            excel_path = tmp_pdf_path.replace(".pdf", ".xlsx")
            with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
                for i, df in enumerate(tables):
                    df.to_excel(writer, sheet_name=f"Sayfa_{i+1}", index=False)

            with open(excel_path, "rb") as f:
                st.success("Dönüştürme tamamlandı! Aşağıdan Excel dosyasını indirebilirsiniz.")
                st.download_button("📥 Excel Dosyasını İndir", f, file_name="donusturulmus_tablo.xlsx")

    except Exception as e:
        st.error(f"Hata oluştu: {str(e)}")

    os.unlink(tmp_pdf_path)
