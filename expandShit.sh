#!/bin/tcsh

set matrix = "parabolic_fem CurlCurl_3 offshore nlpkkt80 CO gsm_106857 msdoor bmw3_2 BenElechi1 t3dh F2 consph SiO2 torso1 dielFilterV3real m_t1 crankseg_2 nd24k mouse_gene human_gene1" 

foreach m ($matrix)
  echo "Expanding $m"
  ./symmetricExpander -i ../sample_matrices/${m}.mtx -o ../sample_matrices/${m}_expanded.mtx
end
