#!/bin/bash
set -e
pegasus_lite_version_major="5"
pegasus_lite_version_minor="0"
pegasus_lite_version_patch="0"
pegasus_lite_enforce_strict_wp_check="true"
pegasus_lite_version_allow_wp_auto_download="true"


. pegasus-lite-common.sh

pegasus_lite_init

# cleanup in case of failures
trap pegasus_lite_signal_int INT
trap pegasus_lite_signal_term TERM
trap pegasus_lite_unexpected_exit EXIT

printf "\n########################[Pegasus Lite] Setting up workdir ########################\n"  1>&2
# work dir
export pegasus_lite_work_dir=$PWD
pegasus_lite_setup_work_dir

printf "\n##############[Pegasus Lite] Figuring out the worker package to use ##############\n"  1>&2
# figure out the worker package to use
pegasus_lite_worker_package

pegasus_lite_section_start stage_in
printf "\n##################### Setting the xbit for executables staged #####################\n"  1>&2
# set the xbit for any executables staged
if [ ! -x cat ]; then
    /bin/chmod +x cat
fi

printf "\n##################### Checking file integrity for input files #####################\n"  1>&2
# do file integrity checks
pegasus-integrity --print-timings --verify=stdin 1>&2 << 'eof'
blastall_00000180_output.txt:blastall_00000077_output.txt:blastall_00000067_output.txt:blastall_00000190_output.txt:blastall_00000170_output.txt:blastall_00000140_output.txt:blastall_00000037_output.txt:blastall_00000130_output.txt:blastall_00000120_output.txt:blastall_00000047_output.txt:blastall_00000057_output.txt:blastall_00000196_output.txt:blastall_00000051_output.txt:blastall_00000041_output.txt:blastall_00000031_output.txt:blastall_00000021_output.txt:blastall_00000150_output.txt:blastall_00000160_output.txt:blastall_00000011_output.txt:blastall_00000097_output.txt:blastall_00000087_output.txt:blastall_00000116_output.txt:blastall_00000106_output.txt:blastall_00000136_output.txt:blastall_00000126_output.txt:blastall_00000186_output.txt:blastall_00000100_output.txt:blastall_00000110_output.txt:blastall_00000107_output.txt:blastall_00000176_output.txt:blastall_00000146_output.txt:blastall_00000137_output.txt:blastall_00000027_output.txt:blastall_00000127_output.txt:blastall_00000166_output.txt:blastall_00000117_output.txt:blastall_00000007_output.txt:blastall_00000017_output.txt:blastall_00000156_output.txt:blastall_00000111_output.txt:blastall_00000101_output.txt:blastall_00000018_output.txt:blastall_00000096_output.txt:blastall_00000008_output.txt:blastall_00000161_output.txt:blastall_00000058_output.txt:blastall_00000078_output.txt:blastall_00000038_output.txt:blastall_00000151_output.txt:blastall_00000141_output.txt:blastall_00000068_output.txt:blastall_00000028_output.txt:blastall_00000131_output.txt:blastall_00000048_output.txt:blastall_00000121_output.txt:blastall_00000098_output.txt:blastall_00000088_output.txt:blastall_00000171_output.txt:blastall_00000191_output.txt:blastall_00000181_output.txt:blastall_00000118_output.txt:blastall_00000125_output.txt:blastall_00000145_output.txt:blastall_00000105_output.txt:blastall_00000138_output.txt:blastall_00000178_output.txt:blastall_00000102_output.txt:blastall_00000009_output.txt:blastall_00000158_output.txt:blastall_00000085_output.txt:blastall_00000185_output.txt:blastall_00000165_output.txt:blastall_00000045_output.txt:blastall_00000065_output.txt:blastall_00000016_output.txt:blastall_00000076_output.txt:blastall_00000036_output.txt:blastall_00000056_output.txt:blastall_00000070_output.txt:blastall_00000187_output.txt:blastall_00000167_output.txt:blastall_00000030_output.txt:blastall_00000147_output.txt:blastall_00000063_output.txt:blastall_00000090_output.txt:blastall_00000010_output.txt:blastall_00000005_output.txt:blastall_00000025_output.txt:blastall_00000129_output.txt:blastall_00000083_output.txt:blastall_00000109_output.txt:blastall_00000050_output.txt:blastall_00000052_output.txt:blastall_00000169_output.txt:blastall_00000149_output.txt:blastall_00000032_output.txt:blastall_00000072_output.txt:blastall_00000092_output.txt:blastall_00000012_output.txt:blastall_00000198_output.txt:blastall_00000061_output.txt:blastall_00000023_output.txt:blastall_00000081_output.txt:blastall_00000189_output.txt:blastall_00000003_output.txt:blastall_00000044_output.txt:blastall_00000024_output.txt:blastall_00000034_output.txt:blastall_00000173_output.txt:blastall_00000004_output.txt:blastall_00000163_output.txt:blastall_00000084_output.txt:blastall_00000094_output.txt:blastall_00000014_output.txt:blastall_00000074_output.txt:blastall_00000183_output.txt:blastall_00000193_output.txt:blastall_00000064_output.txt:blastall_00000054_output.txt:blastall_00000134_output.txt:blastall_00000144_output.txt:blastall_00000104_output.txt:blastall_00000194_output.txt:blastall_00000114_output.txt:blastall_00000124_output.txt:blastall_00000184_output.txt:blastall_00000174_output.txt:blastall_00000154_output.txt:blastall_00000164_output.txt:blastall_00000059_output.txt:blastall_00000039_output.txt:blastall_00000152_output.txt:blastall_00000172_output.txt:blastall_00000069_output.txt:blastall_00000142_output.txt:blastall_00000182_output.txt:blastall_00000103_output.txt:blastall_00000029_output.txt:blastall_00000112_output.txt:blastall_00000192_output.txt:blastall_00000132_output.txt:blastall_00000122_output.txt:blastall_00000049_output.txt:blastall_00000143_output.txt:blastall_00000133_output.txt:blastall_00000153_output.txt:blastall_00000099_output.txt:blastall_00000113_output.txt:blastall_00000079_output.txt:blastall_00000123_output.txt:blastall_00000089_output.txt:blastall_00000162_output.txt:blastall_00000108_output.txt:blastall_00000128_output.txt:blastall_00000115_output.txt:blastall_00000148_output.txt:blastall_00000168_output.txt:blastall_00000135_output.txt:blastall_00000019_output.txt:blastall_00000095_output.txt:blastall_00000006_output.txt:blastall_00000195_output.txt:blastall_00000075_output.txt:blastall_00000175_output.txt:blastall_00000026_output.txt:blastall_00000055_output.txt:blastall_00000155_output.txt:blastall_00000086_output.txt:blastall_00000066_output.txt:blastall_00000046_output.txt:blastall_00000177_output.txt:blastall_00000080_output.txt:blastall_00000157_output.txt:blastall_00000020_output.txt:blastall_00000053_output.txt:blastall_00000015_output.txt:blastall_00000035_output.txt:blastall_00000073_output.txt:cat:blastall_00000040_output.txt:blastall_00000197_output.txt:blastall_00000060_output.txt:blastall_00000093_output.txt:blastall_00000119_output.txt:blastall_00000062_output.txt:blastall_00000179_output.txt:blastall_00000139_output.txt:blastall_00000002_output.txt:blastall_00000022_output.txt:blastall_00000082_output.txt:blastall_00000188_output.txt:blastall_00000159_output.txt:blastall_00000071_output.txt:blastall_00000033_output.txt:blastall_00000091_output.txt:blastall_00000013_output.txt
eof

pegasus_lite_section_end stage_in
set +e
job_ec=0
pegasus_lite_section_start task_execute
printf "\n######################[Pegasus Lite] Executing the user task ######################\n"  1>&2
pegasus-kickstart  -n cat -N cat_00000043 -R condorpool  -s cat_00000043_output.txt=cat_00000043_output.txt -L Blast-Benchmark -T 2021-10-16T00:21:36+00:00 ./cat cat cat_00000043 --path-lock=/tmp/cores.txt.lock --path-cores=/tmp/cores.txt --percent-cpu=0.5 --cpu-work=100 --data --file-size=50 --out=cat_00000043_output.txt
job_ec=$?
pegasus_lite_section_end task_execute
set -e
pegasus_lite_section_start stage_out
pegasus_lite_section_end stage_out

set -e


# clear the trap, and exit cleanly
trap - EXIT
pegasus_lite_final_exit

