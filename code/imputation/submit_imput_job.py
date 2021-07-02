#=========================================
# submit job to TOPMed Imputation Server
#=========================================

import requests
import json

# imputation server url
url = 'https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2'

# add token to header (see Authentication)
headers = {'X-Auth-Token' : 'eyJjdHkiOiJ0ZXh0XC9wbGFpbiIsImFsZyI6IkhTMjU2In0.eyJtYWlsIjoibGluY29sbi5zaGFkZUB1a3kuZWR1IiwiZXhwaXJlIjoxNjI2MzE4ODQ0MDE0LCJuYW1lIjoiTGluY29sbiBTaGFkZSIsImFwaSI6dHJ1ZSwidXNlcm5hbWUiOiJsaW5jb2xuIn0.gSzhyvCkouUuzeWAHlqNMe0_Hbmu6J-C1TEVrwFFOaI' }

# submit new job
vcf = '/home/lmsh224/projects/np_phewas/data/tmp/ADC2_NHW.qced_extract-updated-chr22.vcf.gz';
files = {
  'input-files' : open(vcf, 'rb'),
  'input-refpanel' : 'apps@topmed-r2@1.0.0'
}
r = requests.post(url + '/jobs/submit/imputationserver', files=files, headers=headers)
if r.status_code != 200:
    raise Exception('POST /jobs/submit/imputationserver {}'.format(r.status_code))

# print message
print(r.json()['message'])
print(r.json()['id'])
