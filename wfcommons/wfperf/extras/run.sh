  
#!/bin/bash
for i in {1..9}; do echo "======= $i"; python ten_samples.py frequency -m 14990000 -p 0.$i -s /home/tgcoleman/wfcommons/wfcommons/wfperf/new_test/frequency; done