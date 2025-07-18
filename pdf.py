import streamlit as st
import pdfplumber
import pandas as pd

st.title("PDF'den Tablo Okuma - pdfplumber ile")

uploaded_file = st.file_uploader("PDF dosyanızı yükleyin", type=["pdf"])

if uploaded_file is not None:
    with pdfplumber.open(uploaded_file) as pdf:
        all_tables = []
        for i, page in enumerate(pdf.pages):
            tables = page.extract_tables()
            if tables:
                st.write(f"Sayfa {i+1} tablosu:")
                for table in tables:
                    df = pd.DataFrame(table[1:], columns=table[0])
                    st.dataframe(df)
                    all_tables.append(df)
        if not all_tables:
            st.warning("PDF'de tablo bulunamadı veya PDF taranmış bir görüntü içeriyor.")
