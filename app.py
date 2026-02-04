from flask import Flask, request
from Bio import SeqIO
import os

app = Flask(__name__)

# FASTA 파일 경로 (압축 풀어서 나온 fasta 파일을 지정하세요)
FASTA_FILE = "d1_fasta_clean.fasta"

# 아미노산 분자량 대략값 (단순 평균값, 더 정밀하게 하려면 BioPython의 ProtParam 사용 가능)
aa_weights = {
    'A': 89.1, 'R': 174.2, 'N': 132.1, 'D': 133.1,
    'C': 121.2, 'E': 147.1, 'Q': 146.2, 'G': 75.1,
    'H': 155.2, 'I': 131.2, 'L': 131.2, 'K': 146.2,
    'M': 149.2, 'F': 165.2, 'P': 115.1, 'S': 105.1,
    'T': 119.1, 'W': 204.2, 'Y': 181.2, 'V': 117.1
}

def calc_mw(seq):
    return sum(aa_weights.get(aa, 0) for aa in seq)

@app.route("/search")
def search():
    min_mw = float(request.args.get("min_mw", 10000))
    max_mw = float(request.args.get("max_mw", 50000))
    keyword = request.args.get("keyword", "").lower()

    results = []
    for record in SeqIO.parse(FASTA_FILE, "fasta"):
        mw = calc_mw(str(record.seq))
        if min_mw <= mw <= max_mw:
            if keyword in record.description.lower():
                results.append((record.id, mw, str(record.seq)[:50] + "..."))

    html = "<h2>검색 결과</h2><table border=1><tr><th>ID</th><th>MW</th><th>Sequence (앞 50aa)</th></tr>"
    for rid, mw, seq in results:
        html += f"<tr><td>{rid}</td><td>{mw:.1f}</td><td>{seq}</td></tr>"
    html += "</table>"
    return html

if __name__ == "__main__":
    app.run(debug=True)
