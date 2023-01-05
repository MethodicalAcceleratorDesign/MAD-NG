-- Two dashes start a one-line comment.
--[[
     Adding two ['s and ]'s makes it a
     multi-line comment.
--]]

----------------------------------------------------
-- 1. Variables and flow control.
----------------------------------------------------

[vars]
n = 42  -- All numbers are doubles, but the JIT may specialize them.
-- IEEE-754 64-bit doubles have 52 bits for storing exact int values;
-- machine precision is not a problem for ints < 1e16.

s = 'walternate'  -- Immutable strings like Python.
t = "double-quotes are also fine"
u = [[ Double brackets
       start and end
       multi-line strings.]]
v = "double-quotes \z

     are also fine" -- \z eats next whitespaces
t, u, v = nil  -- Undefines t, u, v.
-- Lua has multiple assignments and nil completion.
-- Lua has garbage collection.

-- Undefined variables return nil. This is not an error:
foo = anUnknownVariable  -- Now foo = nil.
[vars-end]

[cflow]
-- Blocks are denoted with keywords like do/end:
while n < 50 do
  n = n + 1  -- No ++ or += type operators.
end

-- If clauses:
if n > 40 then
  print('over 40')
elseif s ~= 'walternate' then  -- ~= is not equals.
  -- Equality check is == like Python; ok for strs.
  io.write('not over 40\n')  -- Defaults to stdout.
else
  -- Variables are global by default.
  thisIsGlobal = 5  -- Camel case is common.
  -- How to make a variable local:
  local line = io.read()  -- Reads next stdin line.
  -- String concatenation uses the .. operator:
  print('Winter is coming, '..line)
end

-- Only nil and false are falsy; 0 and '' are true!
aBoolValue = false
if not aBoolValue then print('was false') end

-- 'or' and 'and' are short-circuited.
-- This is similar to the a?b:c operator in C/js:
ans = aBoolValue and 'yes' or 'no'  --> ans = 'no'

-- numerical for begin, end[, step] (end included)
revSum = 0
for j = 100, 1, -1 do revSum = revSum + j end
[cflow-end]

----------------------------------------------------
-- 2. Functions.
----------------------------------------------------

[funs]
function fib(n)
  if n < 2 then return 1 end
  return fib(n - 2) + fib(n - 1)
end

-- Closures and anonymous functions are ok:
function adder(x)
  -- The returned function is created when adder is
  -- called, and captures the value of x:
  return function (y) return x + y end
end
a1 = adder(9)
a2 = adder(36)
print(a1(16))  --> 25
print(a2(64))  --> 100

-- Returns, func calls, and assignments all work with lists
-- that may be mismatched in length.
-- Unmatched receivers get nil; unmatched senders are discarded.

x, y, z = 1, 2, 3, 4
-- Now x = 1, y = 2, z = 3, and 4 is thrown away.

function bar(a, b, c)
  print(a, b, c)
  return 4, 8, 15, 16, 23, 42
end

x, y = bar('zaphod')  --> prints "zaphod  nil nil"
-- Now x = 4, y = 8, values 15,..,42 are discarded.

-- Functions are first-class, may be local/global.
-- These are the same:
function f(x) return x * x end
f = function (x) return x * x end

-- And so are these:
local function g(x) return math.sin(x) end
local g; g  = function (x) return math.sin(x) end
-- the 'local g' decl makes g-self-references ok.

-- Calls with one string param don't need parens:
print 'hello'  -- Works fine.
[funs-end]

----------------------------------------------------
-- 3. Tables.
----------------------------------------------------

[tbls]
-- Tables = Lua's only compound data structure;
--   they are associative arrays, i.e. hash-lookup dicts;
--   they can be used as lists, i.e. sequence of non-nil values.

-- Dict literals have string keys by default:
t = {key1 = 'value1', key2 = false, ['key.3'] = true }

-- String keys looking as identifier can use dot notation:
print(t.key1, t['key.3']) -- Prints 'value1 true'.
-- print(t.key.3)         -- Error, needs explicit indexing by string
t.newKey = {}             -- Adds a new key/value pair.
t.key2 = nil              -- Removes key2 from the table.

-- Literal notation for any (non-nil) value as key:
u = {['@!#'] = 'qbert', [{}] = 1729, [6.28] = 'tau'}
print(u[6.28])  -- prints "tau"

-- Key matching is basically by value for numbers
-- and strings, but by identity for tables.
a = u['@!#']  -- Now a = 'qbert'.
b = u[{}]     -- We might expect 1729, but it's nil:

-- A one-table-param function call needs no parens:
function h(x) print(x.key1) end
h{key1 = 'Sonmi~451'}  -- Prints 'Sonmi~451'.

for key, val in pairs(u) do -- Table iteration.
  print(key, val)
end

-- List literals implicitly set up int keys:
l = {'value1', 'value2', 1.21, 'gigawatts'}
for i,v in ipairs(l) do  -- List iteration.
  print(i,v,l[i])        -- Indices start at 1 !
end
print("length=", #l)     -- # is defined only for sequence.
-- A 'list' is not a real type, l is just a table
-- with consecutive integer keys, treated as a list,
-- i.e. l = {[1]='value1', [2]='value2', [3]=1.21, [4]='gigawatts'}
-- A 'sequence' is a list with non-nil values.
[tbls-end]

----------------------------------------------------
-- 4. Methods.
----------------------------------------------------

[mthds]
-- Methods notation:
--   function tblname:fn(...) is the same as
--     function tblname.fn(self, ...) with self being the table.
--   calling tblname:fn(...) is the same as
--     tblname.fn(tblname, ...)       here self becomes the table.
t = { disp=function(s) print(s.msg) end, -- Method 'disp'
      msg="Hello world!" }
t:disp() -- Prints "Hello world!"
function t:setmsg(msg) self.msg=msg end  -- Add a new method 'setmsg'
t:setmsg "Good bye!"
t:disp() -- Prints "Good bye!"
[mthds-end]
