# FORAL

## Fingerprint

Bi programa daude:
  - **Fingerprint1.f90**: Honek masa molekularra biderkatzen dio atomo bakoitzari *d* dimentsioko fingerprint-a sortuz.
  - **Fingerprint2.f90**: Honek *2d* dimenstioko fingerprint-a sortzen du. Lehenengo zatia *Si* atomoak bakarrik kontuan hartuta eta bigarren zatia *O* atomoak bakarrik kontuan hartuta.
  
  *Gogoratu kasu bietan fingerprint-ari balio bat gehitu zaiola zero puntuan atomo horren masa molekularrekin*

### Konpilatzeko:

Programa biak *gfortran*-ekin konpilatzen dira hurrengo komandoarekin:

> gfortran fingerprint1.f90 -o finger1

### Exekutatzeko:

> ./finger1 [atomo_mota_kop] [d] [R_cut] [input_file] [output_file]

Edozein duda eduki ezekero kontaktatu nirekin *xabier.mendez@ehu.eus* helbide elektronikoan.
