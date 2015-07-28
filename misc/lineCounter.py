f=open("proband.txt", 'r+')
counter = 0
for line in f:
	counter +=1
print str(counter) + " paired mates skipped."