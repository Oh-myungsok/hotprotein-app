import streamlit as st
from Bio import SeqIO

FASTA_FILE = "d1_fasta_clean.fasta"

st.title("Protein Calculator & Search")

# 입력 영역
min_mw = st.number_input("Minimum molecular weight", value=10000)
max_mw = st.number_input("Maximum molecular weight", value=50000)
keyword = st.text_input("Keyword", value="DNA")

# 버튼 클릭 시 검색 실행
if st.button("Search"):
    results = []
    for record in SeqIO.parse(FASTA_FILE, "fasta"):
        desc = record.description
        # 단순 분자량 계산 (아미노산 평균 110 Da 가정)
        mw = len(record.seq) * 110
        if min_mw <= mw <= max_mw and keyword.lower() in desc.lower():
            results.append({"ID": record.id, "MW": mw, "Description": desc})

    if results:
        st.write("### Search Results")
        st.dataframe(results)
    else:
        st.warning("No matches found.")
