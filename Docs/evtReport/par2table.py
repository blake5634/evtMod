import sys

for line in sys.stdin:
   line = line.replace('_','\_')
   star = False
   if '*' in line:
      star = True
      line = line.replace('*','')
   pts = line.split()
   par = pts[0]
   if star:
      par += '*'
   val = pts[1]
   units = pts[2:]
   units = ''.join(units)
   print(f'{par}  &  {val} & {units} \\\\', file=sys.stdout)

