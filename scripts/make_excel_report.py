import re
import os
import gzip
from sys import argv
from typing import Dict, List

import xlsxwriter
from Bio import SeqIO


def parser_result_files(file: str) -> Dict[str, dict]:
    d = {'results': {}, 'params': {}}
    with open(file, 'r') as infile:
        for line in infile:
            prefix, key, value = line.strip().split("\t")
            d[prefix][key] = value
    return d


def base_info(d: Dict[str, dict], book: xlsxwriter.workbook.Workbook, row=1, ont=True):
    sheet1 = book.get_worksheet_by_name('基本信息')
    sheet2 = book.get_worksheet_by_name('统计图表')
    bold_field = book.add_format()
    bold_field.set_bold()
    params = d['params']
    results = d['results']
    sheet1.set_column('A:A', 10)
    sheet1.set_column('B:B', 35)
    sheet1.set_column('C:E', 26)
    sheet1.write(f"A{row}", '基本信息')
    row += 1
    for para_key, para_value in params.items():
        sheet1.write_string(f"B{row}", para_key, bold_field)
        sheet1.write_string(f"C{row}", para_value)
        row += 1
    sheet1.write(f"A{row}", "数据信息")
    row += 1
    if ont:
        with open(results['summary']) as infile:
            first_line = next(infile)
            _, analysis_name, sample_name = first_line.strip().split("\t")
            sheet1.write_string(f"B{row}", "analysis_name", bold_field)
            sheet1.write_string(f"C{row}", analysis_name)
            row += 1
            sheet1.write_string(f"B{row}", "sample_name", bold_field)
            sheet1.write_string(f"C{row}", sample_name)
            row += 1
            for line in infile:
                field, raw, clean, mapped = line.strip().split("\t")
                sheet1.write_string(f"B{row}", field, bold_field)
                sheet1.write_string(f"C{row}", raw)
                sheet1.write_string(f"D{row}", clean)
                sheet1.write_string(f"E{row}", mapped)
                row += 1
        for new_row, title, png in zip([1, 28, 55],
                                       ['RawData', 'CleanData', 'MappedData'],
                                       [results['rawPng'], results['cleanPng'], results['mappedPng']]):
            sheet2.write(f"A{new_row}", f"统计图表-{title}", bold_field)
            sheet2.insert_image(f"B{new_row + 1}", png)

        sheet2.set_column('M:M', 12)
        for new_row, title, png in zip([1, 30],
                                       ['原始深度统计', '1000X深度统计'],
                                       [results['rawCovPng'], results['modifiedCovPng']]):
            sheet2.write(f"M{new_row}", title, bold_field)
            sheet2.insert_image(f"N{new_row + 1}", png, {'x_scale': 0.9, 'y_scale': 0.9})

    else:
        # For ngs excel report base info
        with open(results['summary']) as infile:
            for line in infile:
                field, number = line.strip().split("\t")
                sheet1.write_string(f"B{row}", field, bold_field)
                sheet1.write_string(f"C{row}", number)
                row += 1
        sheet2.set_column('A:A', 12)
        for new_row, title, png in zip([1, 30],
                                       ['原始深度统计', '1000X深度统计'],
                                       [results['rawCovPng'], results['modifiedCovPng']]):
            sheet2.write(f"A{new_row}", title, bold_field)
            sheet2.insert_image(f"B{new_row + 1}", png, {'x_scale': 0.9, 'y_scale': 0.9})


def pango2excel(pango: str, book: xlsxwriter.workbook.Workbook, row=1) -> int:
    if os.path.getsize(pango) == 0:
        return row
    bold_field = book.add_format()
    bold_field.set_bold()
    sheet = book.get_worksheet_by_name('序列及变异')
    sheet.set_column("A:A", 15)
    sheet.write_string(f"A{row}", "Pangolin结果", bold_field)
    row += 1
    with open(pango, 'r') as infile:
        for line in infile:
            col = 1
            for field in line.strip().split(","):
                sheet.write_string(row - 1, col, field)
                col += 1
            row += 1
    return row


def write_consensus(consensus, book: xlsxwriter.workbook.Workbook, row) -> int:
    sheet = book.get_worksheet_by_name('序列及变异')
    bold_field = book.add_format()
    bold_field.set_bold()
    sheet.write_string(f"A{row}", '一致性序列', bold_field)
    row += 1
    for record in SeqIO.parse(open(consensus, 'r'), 'fasta'):
        sheet.write_string(f"B{row}", record.id)
        row += 1
        sheet.write_string(f"B{row}", str(record.seq))
        row += 1
    return row


def vcf2list(vcf: str, ont=True) -> List[list]:
    snpeff = False
    header_with_snpeff = ['染色体', '位置', '参考碱基', '变异碱基', '深度', 'Ref深度', 'Alt深度', '基因', '核酸变化',
                          '氨基酸变化']
    header_without_snpeff = ['染色体', '位置', '参考碱基', '变异碱基', '深度', 'Ref深度', 'Alt深度']
    header = header_without_snpeff
    rows_content = []
    with gzip.open(vcf, 'rt') if vcf.endswith('vcf.gz') else open(vcf, 'r') as infile:
        for line in infile:
            if line.startswith("#"):
                if line.startswith("##SnpEffCmd"):
                    snpeff = True
                if line.startswith("##INFO=<ID=ANN,") and snpeff:
                    try:
                        mat = re.search("\"Functional annotations: '(.+)'.+\"", line)
                        snpeff_columns = mat.group(1).split(" | ")
                        gene_name_idx = snpeff_columns.index('Gene_Name')
                        hgvs_c_idx = snpeff_columns.index('HGVS.c')
                        hgvs_p_idx = snpeff_columns.index('HGVS.p')
                        header = header_with_snpeff
                    except (AttributeError, ValueError) as error:
                        header = header_without_snpeff
                continue
            contig, pos, _, ref, alt, _, _, info, fmt, sample = line.strip('\n').split("\t")
            sample_dict = {x: y for x, y in zip(fmt.split(":"), sample.split(":"))}

            info_dict = {k: v for k, v in [each_info.split('=') for each_info in info.split(';') if "=" in each_info]}

            dp = sample_dict['DP']
            ref_dp, alt_dp = sample_dict['AD'].split(',')
            row_content = [contig, pos, ref, alt, dp, ref_dp, alt_dp]
            if header == header_with_snpeff:
                ann = info_dict['ANN']
                ann_fields = ann.split('|')
                gene_name = ann_fields[gene_name_idx]
                hgvs_c = ann_fields[hgvs_c_idx]
                hgvs_p = ann_fields[hgvs_p_idx]
                row_content = [contig, pos, ref, alt, dp, ref_dp, alt_dp, gene_name, hgvs_c, hgvs_p]
                row_content = list(map(lambda x: x if x else '-', row_content))
            rows_content.append(row_content)
    rows_content.insert(0, header)
    return rows_content


def vcf2excel(var_list: List[list], book: xlsxwriter.workbook.Workbook, row: int):
    # write excel with rows_content
    sheet = book.get_worksheet_by_name('序列及变异')
    bold_field = book.add_format()
    bold_field.set_bold()
    sheet.write_string(f"A{row}", '变异信息', bold_field)
    for var in var_list:
        col = 1
        for field in var:
            sheet.write_string(row, col, field)
            col += 1
        row += 1
    return 20


def main():
    infile = argv[1]
    out_excel = argv[2]
    d = parser_result_files(infile)
    ont = d['params']['long_reads'] != 'null'
    book = xlsxwriter.Workbook(out_excel)
    sheet1 = book.add_worksheet('基本信息')
    sheet2 = book.add_worksheet('统计图表')
    sheet3 = book.add_worksheet('序列及变异')
    base_info(d, book, ont=ont)
    row = pango2excel(d['results']['pangoRpt'], book)
    row = write_consensus(d['results']['consensus'], book, row=row)
    var_list = vcf2list(d['results']['vcfRpt'], ont=ont)
    vcf2excel(var_list, book, row=row)
    book.close()


if __name__ == '__main__':
    main()
