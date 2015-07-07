nplusone = 8

def getNext(n):
	s = ('{0:0' + str(nplusone - 1) + 'b}').format(n);
	s = list(s)
	s = map(lambda x: int(x), s)
	return s

def findIt():
	for i in range(0, 2**nplusone):
		if(i % 25000000 == 0):
			print i
		l = getNext(i)

		if sum(l) != nplusone/2 and sum(l) != nplusone/2 + 1:
			continue

		l.extend(l)

		tracker = [0] * nplusone
		tracker1 = [0] * nplusone
		for start in range(0, nplusone - 1):
			for cur in range(1, nplusone - 1):
				if l[start + cur] == 1:
					tracker[cur] = tracker[cur] + 1
					tracker1[cur] = 1
		if sum(tracker) == nplusone - 2 and sum(tracker1) == nplusone-2:
			return l;

	return "could not find it"

print findIt()[:nplusone]