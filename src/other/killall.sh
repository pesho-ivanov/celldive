kill -9 `ps -ef | grep bowtie | grep -v grep | awk '{print $2}'`
kill -9 `ps -ef | grep kallisto | grep -v grep | awk '{print $2}'`
