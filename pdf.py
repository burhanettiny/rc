import streamlit as st
import pdfplumber
import pandas as pd
from io import BytesIO
from PIL import Image
import pytesseract

# Eğer Tesseract yüklü değilse, path ayarlanmalı
# pytesseract.pytesseract.tesseract_cmd = r"C:\Program Files\Tesseract-OCR\tesseract.exe"

st.title("PDF'den Excel'e (Dijital ve El Yazısı Tablolar)")

uploaded_file = st.file_uploader("PDF dosyanızı yükleyin", type=["pdf"])

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

def ocr_to_dataframe(image):
    # OCR verisini dataframe olarak al
    df_data = pytesseract.image_to_data(image, output_type=pytesseract.Output.DATAFRAME)
    df_data = df_data[df_data.conf != -1]  # Geçerli kelimeler
    
    if df_data.empty:
        return None
    
    # Satır gruplama
    rows = []
    for top_val, group in df_data.groupby('top'):
        row = list(group.sort_values('left')['text'])
        rows.append(row)
    
    if len(rows) > 1:
        try:
            df = pd.DataFrame(rows[1:], columns=make_unique_columns(rows[0]))
        except:
            df = pd.DataFrame(rows)
    else:
        df = pd.DataFrame(rows)
    
    return df

if uploaded_file is not None:
    all_tables = []
    with pdfplumber.open(uploaded_file) as pdf:
        for i, page in enumerate(pdf.pages):
            # 1️⃣ Dijital tablolardan alma
            tables = page.extract_tables()
            if tables:
                for j, table in enumerate(tables):
                    df = pd.DataFrame(table[1:], columns=table[0])
                    df.columns = make_unique_columns(df.columns)
                    all_tables.append((f"Sayfa {i+1} Tablo {j+1} (Dijital)", df))
            
            # 2️⃣ OCR ile el yazısı/görsel tablolar
            im = page.to_image(resolution=300).original
            df_ocr = ocr_to_dataframe(im)
            if df_ocr is not None and not df_ocr.empty:
                all_tables.append((f"Sayfa {i+1} Tablo (OCR)", df_ocr))
    
    if all_tables:
        export_option = st.radio("Excel dosyasını nasıl oluşturmak istersiniz?", 
                                 ("Her tablo ayrı sayfada", "Tüm tablolar tek sayfada"))

        output = BytesIO()
        with pd.ExcelWriter(output, engine="xlsxwriter") as writer:
            if export_option == "Her tablo ayrı sayfada":
                for sheet_name, df in all_tables:
                    safe_sheet_name = sheet_name[:31]
                    df.to_excel(writer, sheet_name=safe_sheet_name, index=False)
            else:
                worksheet_name = "Tüm Tablolar"
                worksheet = writer.book.add_worksheet(worksheet_name)
                writer.sheets[worksheet_name] = worksheet
                startrow = 0
                for title, df in all_tables:
                    worksheet.write(startrow, 0, title)
                    df.to_excel(writer, sheet_name=worksheet_name, startrow=startrow + 1, index=False, header=True)
                    startrow += len(df) + 3

        processed_data = output.getvalue()

        st.download_button(
            label="Excel Dosyasını İndir",
            data=processed_data,
            file_name="pdf_tablosu.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )
    else:
        st.warning("PDF'de tablo bulunamadı veya OCR ile okunabilecek tablo yok.")
