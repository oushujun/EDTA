hmmpress -h

for hmm in database/*hmm
do
	hmmpress -f $hmm
done
