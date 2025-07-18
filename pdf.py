import streamlit as st
import pdfplumber
import pandas as pd
from io import BytesIO

def make_unique_columns(columns):
    seen = {}
    result = []
    for col in columns:
        if col not in seen:
            seen[col] = 0
            result.append(col)
        else:
            seen[col] += 1
            result.append(f"{col}_{seen[col]}")
    return result

st.title("PDF' to excel")

uploaded_file = st.file_uploader("PDF dosyanızı yükleyin", type=["pdf"])

if uploaded_file is not None:
    all_tables = []
    with pdfplumber.open(uploaded_file) as pdf:
        for i, page in enumerate(pdf.pages):
            tables = page.extract_tables()
            if tables:
                st.write(f"Sayfa {i+1} tablosu:")
                for j, table in enumerate(tables):
                    df = pd.DataFrame(table[1:], columns=table[0])
                    df.columns = make_unique_columns(df.columns)
                    st.dataframe(df)
                    all_tables.append((f"Sayfa_{i+1}_Tablo_{j+1}", df))

    if all_tables:
        # Excel dosyasını oluştur
        output = BytesIO()
        with pd.ExcelWriter(output, engine="xlsxwriter") as writer:
            for sheet_name, df in all_tables:
                # Excel sayfa adı 31 karakterden uzun olmamalı, onu kontrol edelim
                safe_sheet_name = sheet_name[:31]
                df.to_excel(writer, sheet_name=safe_sheet_name, index=False)
            with pd.ExcelWriter(output, engine="xlsxwriter") as writer:
                for sheet_name, df in all_tables:
                    safe_sheet_name = sheet_name[:31]
                    df.to_excel(writer, sheet_name=safe_sheet_name, index=False)
            processed_data = output.getvalue()

        # Excel dosyasını oluştur
        output = BytesIO()
        with pd.ExcelWriter(output, engine="xlsxwriter") as writer:
            for sheet_name, df in all_tables:
                # Excel sayfa adı 31 karakterden uzun olmamalı, onu kontrol edelim
                safe_sheet_name = sheet_name[:31]
                df.to_excel(writer, sheet_name=safe_sheet_name, index=False)
            with pd.ExcelWriter(output, engine="xlsxwriter") as writer:
                for sheet_name, df in all_tables:
                    safe_sheet_name = sheet_name[:31]
                    df.to_excel(writer, sheet_name=safe_sheet_name, index=False)
            processed_data = output.getvalue()    

        st.download_button(
            label="Excel Dosyasını İndir",
            data=processed_data,
            file_name="pdf_tablosu.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )
    else:
        st.warning("PDF'de tablo bulunamadı veya PDF taranmış bir görüntü içeriyor.")
