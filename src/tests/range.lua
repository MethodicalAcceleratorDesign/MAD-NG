local jitv   = require('jit.v')
-- jit.opt.start("sizemcode=256","loopunroll=25")
-- jitv.start()

local range = require 'range'
local gmath = require 'gmath'

local isrange = gmath.is_range

local function check(start, stop, step)
  local r = range(start, stop, step)
  print('range=', r:range())
  print('size=', r:size())
  for i,v in ipairs(r) do
    print('i=',i,', v=', v)
  end
  print()
end

check(0)
check(1)
check(-1)

check( 0, 0)
check( 0, 1)
check( 1, 0)

check(-1, 1)
check( 1,-1)
check(-1, 1,-1)
check( 1,-1, 1)

local r=range(1,1.7,0.1)
print('range=', r:range())
for i,v in ipairs(r) do
  print(i,v)
end

local t = r:totable()
print('table=', t)
for i,v in ipairs(t) do
  print(i,v)
end

local x = r:tovector()
print('vector=', x)
for i,v in ipairs(x) do
  print(i,v)
end

print("lengths=", #r, #t, #x)
print("values=", r[0], t[0], x[0])
print("values=", r[1], t[1], x[1])
print("values=", r[#r+1], t[#t+1], x[#x+1])
io.write('r=', string.format("%.14f", r[#r]), '\n')
io.write('t=', string.format("%.14f", t[#t]), '\n')
io.write('x=', string.format("%.14f", x[#x]), '\n')


print(isrange(nil))
print(isrange(2))
print(isrange(range))
print(isrange({}))
print(isrange(r))
