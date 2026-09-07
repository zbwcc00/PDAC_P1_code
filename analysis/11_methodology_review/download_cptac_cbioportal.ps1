$root = 'D:\第二篇大论文'
$out = Join-Path $root 'data\05_protein\CPTAC_PAAD_cBioPortal'
New-Item -ItemType Directory -Force -Path $out | Out-Null
Write-Output "output_dir=$out exists=$(Test-Path $out)"
$study = 'paad_cptac_2021'
$samples = Invoke-RestMethod "https://www.cbioportal.org/api/studies/$study/samples?projection=SUMMARY"
$ids = @($samples.sampleId)
$body = @{entrezGeneIds=@(2034,6867,4082,9709,5175,7450,51705);sampleIds=$ids} | ConvertTo-Json -Depth 4
$protein = Invoke-RestMethod -Method Post -Uri "https://www.cbioportal.org/api/molecular-profiles/${study}_protein_quantification/molecular-data/fetch" -Body $body -ContentType 'application/json'
$protein | Export-Csv (Join-Path $out 'protein_quantification.tsv') -NoTypeInformation -Delimiter "`t"
$clinical = Invoke-RestMethod "https://www.cbioportal.org/api/studies/$study/clinical-data?clinicalDataType=PATIENT&projection=DETAILED"
$clinical | Export-Csv (Join-Path $out 'patient_clinical_data.tsv') -NoTypeInformation -Delimiter "`t"
$samples | Export-Csv (Join-Path $out 'sample_metadata.tsv') -NoTypeInformation -Delimiter "`t"
Write-Output "samples=$($samples.Count); protein_records=$($protein.Count); clinical_records=$($clinical.Count)"
