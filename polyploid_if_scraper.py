# @Version   :1.0
# @Author    :杜伟轩
# @File      :polyploid_if_scraper.py
# @Time      :2025/8/8 12:52
# -*- coding: utf-8 -*-

import argparse
import time
import requests
import csv
from datetime import datetime, timedelta
import pandas as pd
from xml.etree import ElementTree as ET

# -----------------------
# 配置
# -----------------------
CROSSREF_BASE = "https://api.crossref.org/works"
ENTREZ_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
USER_EMAIL = "dwx741886@163.com"  # 替换为你的邮箱
TOOL_NAME = "polyploid_scraper"

# 搜索关键词
kw = []
while True:
    print("输入完毕请填写“None”")
    i = input("请输入关键词：")
    kw.append(i)
    if i == "None":
        break

SEARCH_TERMS = kw
# SEARCH_TERMS = [
#     '"plant polyploidy"',
#     '"plant polyploidization"',
#     '"plant whole genome duplication"',
#     '"crop polyploidy"',
#     '"crop polyploidization"',
#     '"crop whole genome duplication"',
#     '"vegetable polyploidy"',
#     '"vegetable polyploidization"',
#     '"vegetable whole genome duplication"',
#     '"植物 多倍化"',
#     '"植物 多倍体"',
#     '"作物 多倍化"',
#     '"作物 多倍体"',
#     '"蔬菜 多倍化"',
#     '"蔬菜 多倍体"'
# ]

YEARS_BACK = 5
REQUESTS_DELAY = 1.0  # 秒


# -----------------------
# 工具函数
# -----------------------
def five_years_ago_iso():
    dt = datetime.utcnow() - timedelta(days=YEARS_BACK * 365)
    return dt.strftime("%Y-%m-%d")


def normalize_journal_name(name: str) -> str:
    if not name:
        return ""
    return " ".join(name.lower().strip().split())


def load_if_csv(path):
    mapping = {}
    with open(path, newline='', encoding='utf-8') as f:
        reader = csv.reader(f)
        for row in reader:
            if not row: continue
            try:
                if_val = float(row[1])
            except:
                continue
            mapping[normalize_journal_name(row[0])] = if_val
    return mapping


# -----------------------
# 自动分类函数
# -----------------------
def classify_polyploid(keywords, title):
    text = (keywords + " " + title).lower()
    if "allopolyploid" in text or "allopolyploidy" in text or "allopolyploidization" in text:
        return "allopolyploid"
    elif "autopolyploid" in text or "autopolyploidy" in text or "autopolyploidization" in text:
        return "autopolyploid"
    else:
        return "unknown"


# -----------------------
# Crossref 搜索
# -----------------------
def search_crossref(query, rows=100):
    params = {
        "query.bibliographic": query,
        "rows": rows,
        "filter": f"from-pub-date:{five_years_ago_iso()}"
    }
    headers = {"User-Agent": f"{TOOL_NAME} (mailto:{USER_EMAIL})"}
    resp = requests.get(CROSSREF_BASE, params=params, headers=headers, timeout=30)
    time.sleep(REQUESTS_DELAY)
    resp.raise_for_status()
    data = resp.json()
    return data.get("message", {}).get("items", [])


def extract_from_crossref_item(item):
    title = item.get("title", [""])[0] if item.get("title") else ""
    doi = item.get("DOI", "")
    url = item.get("URL", "")
    journal = ""
    if "container-title" in item and item["container-title"]:
        journal = item["container-title"][0]
    pub_date = ""
    if item.get("published-print"):
        parts = item["published-print"].get("date-parts", [[]])
        if parts and parts[0]:
            pub_date = "-".join(str(x) for x in parts[0])
    elif item.get("published-online"):
        parts = item["published-online"].get("date-parts", [[]])
        if parts and parts[0]:
            pub_date = "-".join(str(x) for x in parts[0])
    keywords = item.get("subject", []) or []
    return {
        "title": title,
        "doi": doi,
        "url": url,
        "journal": journal,
        "pub_date": pub_date,
        "keywords": "; ".join(keywords)
    }


# -----------------------
# PubMed 检索
# -----------------------
def search_pubmed(term, retmax=200):
    params = {
        "db": "pubmed",
        "term": term,
        "retmax": retmax,
        "retmode": "xml",
        "datetype": "pdat",
        "mindate": five_years_ago_iso(),
        "tool": TOOL_NAME,
        "email": USER_EMAIL
    }
    resp = requests.get(f"{ENTREZ_BASE}/esearch.fcgi", params=params, timeout=30)
    time.sleep(REQUESTS_DELAY)
    resp.raise_for_status()
    root = ET.fromstring(resp.text)
    return [elem.text for elem in root.findall(".//Id")]


def fetch_pubmed_details(id_list):
    if not id_list:
        return []
    ids = ",".join(id_list)
    params = {
        "db": "pubmed",
        "id": ids,
        "retmode": "xml",
        "tool": TOOL_NAME,
        "email": USER_EMAIL
    }
    resp = requests.get(f"{ENTREZ_BASE}/efetch.fcgi", params=params, timeout=30)
    time.sleep(REQUESTS_DELAY)
    resp.raise_for_status()
    root = ET.fromstring(resp.text)
    results = []
    for article in root.findall(".//PubmedArticle"):
        try:
            art_info = article.find(".//Article")
            title = art_info.findtext("ArticleTitle") or ""
            journal = art_info.findtext(".//Journal/Title") or ""
            pubdate_elem = art_info.find(".//Journal/JournalIssue/PubDate")
            pub_date = ""
            if pubdate_elem is not None:
                y = pubdate_elem.findtext("Year")
                m = pubdate_elem.findtext("Month")
                d = pubdate_elem.findtext("Day")
                parts = [p for p in (y, m, d) if p]
                pub_date = "-".join(parts)
            doi = ""
            for eid in article.findall(".//ArticleId"):
                if eid.get("IdType") == "doi":
                    doi = eid.text
            mesh_terms = [m.findtext("DescriptorName") for m in article.findall(".//MeshHeading") if
                          m.findtext("DescriptorName")]
            url = f"https://doi.org/{doi}" if doi else ""
            results.append({
                "title": title,
                "journal": journal,
                "pub_date": pub_date,
                "doi": doi,
                "keywords": "; ".join(mesh_terms),
                "url": url
            })
        except:
            continue
    return results


# -----------------------
# 主流程
# -----------------------
def main(args):
    if_map = {}
    if args.if_csv:
        print(f"[+] 载入影响因子 CSV：{args.if_csv}")
        if_map = load_if_csv(args.if_csv)
        print(f"    已载入 {len(if_map)} 条记录")

    results = {}
    # Crossref
    for term in SEARCH_TERMS:
        items = search_crossref(term, rows=200)
        for it in items:
            rec = extract_from_crossref_item(it)
            key = rec.get("doi") or rec.get("title")[:200]
            if key not in results:
                results[key] = rec

    # PubMed
    pmids = set()
    for term in SEARCH_TERMS:
        for p in search_pubmed(term, retmax=200):
            pmids.add(p)
    for i in range(0, len(pmids), 100):
        for d in fetch_pubmed_details(list(pmids)[i:i + 100]):
            key = d.get("doi") or d.get("title")[:200]
            if key not in results:
                results[key] = d
            else:
                if d.get("keywords"):
                    results[key]["keywords"] = d.get("keywords")

    # 筛选 IF≥7
    final_list = []
    for rec in results.values():
        jname = normalize_journal_name(rec.get("journal", ""))
        if not jname:
            continue
        if jname in if_map and if_map[jname] >= 7.0:
            rec["journal_if"] = if_map[jname]
            rec["category"] = classify_polyploid(rec.get("keywords", ""), rec.get("title", ""))
            final_list.append(rec)

    # 输出 Excel
    df = pd.DataFrame(final_list)
    df = df[["title", "pub_date", "keywords", "category", "doi", "url", "journal", "journal_if"]]
    df.to_excel(args.out, index=False)
    print(f"[+] 已保存 {len(df)} 条记录到 {args.out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="检索近5年植物多倍化相关论文(IF>=7)")
    parser.add_argument("--out", required=True, help="输出 Excel 文件名")
    parser.add_argument("--if_csv", required=True, help="期刊影响因子 CSV 文件")
    args = parser.parse_args()
    main(args)
