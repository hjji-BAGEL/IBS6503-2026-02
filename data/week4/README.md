# 4주차 실습 데이터

2주차에 찾은 Himes et al. (2014) 덱사메타손 실험의 실제 read 일부입니다. 사람 1번 염색체의 약 330 kb 구간(hg19 chr1:11,053,773–11,386,132)에 붙은 read만 남겼습니다.

| 파일 | 내용 |
|---|---|
| `1_control.R1.fastq.gz`, `1_control.R2.fastq.gz` 외 3샘플 | read 서열과 품질. 63 bp paired-end. R1과 R2는 같은 조각의 양 끝 |
| `hg19_chr1_11.4Mb.fa.gz` | 참조 서열. hg19 chr1의 앞 11,400,000 bp. 앞 11,000,000 bp는 `N`으로 채움 |
| `genes_GRCh37.75_chr1_11Mb.gtf` | 이 구간의 유전자 위치 (Ensembl GRCh37 release 75) |
| `samples.tsv` | 샘플 이름과 원 실험의 run, 세포주, 처리 여부 |

샘플 이름은 `data/airway_scaledcounts.subset.tsv`의 열 이름과 같습니다.

참조 서열의 앞부분을 `N`으로 채운 것은 정렬 결과가 실제 hg19 좌표로 나오게 하려는 것입니다. 그래서 IGV의 hg19 유전체에 그대로 얹을 수 있습니다.

원본은 고치지 않습니다. 정렬 결과는 `results/`에 새로 만듭니다.

## 출처

- Bioconductor `airway` 패키지 1.32.0의 `inst/extdata/*_subset.bam`에서 read를 다시 꺼낸 것. doi:10.18129/B9.bioc.airway
- 원 실험: Himes BE, et al. *RNA-Seq transcriptome profiling identifies CRISPLD2 as a glucocorticoid responsive gene that modulates cytokine function in airway smooth muscle cells.* PLoS One 9(6):e99625 (2014). PMID 24926665. GEO GSE52778
- 참조 서열: UCSC Genome Browser hg19 (GRCh37)
