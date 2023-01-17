import sys, csv


lightWave = []
power = []
with open(sys.argv[1]) as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',')
    firstRow = True
    for row in spamreader:
        if firstRow:
            firstRow = False
            lightWave = row
        else:
            power = row
            break

output = [-1] * (701 - 350)
for i in range(len(lightWave)):
    lw = int(lightWave[i])
    pow = float(power[i])
    output[lw - 350] = pow

# On suppose que la première valeur est bien définie
for i in range(1, len(output)):
    if output[i] >= 0:
        continue
    j = i-1
    while output[j] == -1:
        j=j-1
    k = i+1
    while output[k] == -1:
        k=k+1
    output[i] = output[j] + (i - j) * (output[k] - output[j]) / (k - j)
#print(output)

with open(sys.argv[2], "w") as outputfile:
    for i in range(0, len(output), 20):
        print(i + 350)
        outputfile.write("{}\n".format(output[i]))

