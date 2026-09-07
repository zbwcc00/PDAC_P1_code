$ErrorActionPreference = 'Stop'
$root = 'D:\第二篇大论文'
$input = Join-Path $root 'analysis\08_program_validation\endothelial_target_prioritization.tsv'
$out = Join-Path $root 'analysis\10_communication_pseudotime\Target_screen\HPA_target_summary.tsv'
$mapping = Import-Csv (Join-Path $root 'analysis\10_communication_pseudotime\Target_screen\endothelial_targets_ensembl.tsv') -Delimiter "`t"
$rows = @()
foreach ($item in $mapping) {
  $gene = $item.SYMBOL
  try {
    $ensg = $item.ENSEMBL
    $hpa = Invoke-RestMethod -Uri ("https://www.proteinatlas.org/{0}.json" -f $ensg) -TimeoutSec 20
    $enrich = if ($hpa.'RNA tissue cell type enrichment') { ($hpa.'RNA tissue cell type enrichment' -join ';') } else { '' }
    $proteinClass = if ($hpa.'Protein class') { ($hpa.'Protein class' -join ';') } else { '' }
    $rows += [pscustomobject]@{gene=$gene; ensembl=$ensg; evidence=$hpa.Evidence; hpa_evidence=$hpa.'HPA evidence'; rna_sc_specificity=$hpa.'RNA single cell type specificity'; rna_sc_distribution=$hpa.'RNA single cell type distribution'; protein_sc_specificity=$hpa.'Protein cell type specificity'; protein_sc_distribution=$hpa.'Protein cell type distribution'; tissue_cell_type_enrichment=$enrich; protein_class=$proteinClass; cancer_specificity=$hpa.'RNA cancer specificity'}
  } catch {
    $rows += [pscustomobject]@{gene=$gene; ensembl=''; evidence='query_failed'; hpa_evidence=$_.Exception.Message; rna_sc_specificity=''; rna_sc_distribution=''; protein_sc_specificity=''; protein_sc_distribution=''; tissue_cell_type_enrichment=''; protein_class=''; cancer_specificity=''}
  }
  Start-Sleep -Milliseconds 200
}
$rows | Export-Csv -Path $out -Delimiter "`t" -NoTypeInformation -Encoding UTF8
Write-Output ("Wrote {0}" -f $out)
