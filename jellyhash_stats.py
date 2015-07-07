def median(lst):
    lst = sorted(lst)
    if len(lst) < 1:
            return None
    if len(lst) %2 == 1:
            return lst[((len(lst)+1)/2)-1]
    else:
            return float(sum(lst[(len(lst)/2)-1:(len(lst)/2)+1]))/2.0

infile = open('jellyhash_bucket_histo')

i=0;

bucket_sizes = []
zeros = 0

for line in infile:
	n = int(line.split(" ")[1])
	if i%25000000 == 0:
		print i, n
	if n == 0:
		zeros = zeros + 1
	else:
		bucket_sizes.append(n)
	i += 1

infile.close()

print "length including 0s:", len(bucket_sizes) + zeros
print "number non-zeros:", len(bucket_sizes)
print "NOTE: for below stats, we don't count zeros"
print "min:", min(bucket_sizes)
print "max:", max(bucket_sizes)
print "mean:", sum(bucket_sizes) * 1.0 / len(bucket_sizes)
print "median:", median(bucket_sizes)