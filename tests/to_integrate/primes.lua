local isPrimeDivisible = function(self, c)
  for i=2, self.prime_count do
    if self.primes[i] * self.primes[i] > c then break end
    if c % self.primes[i] == 0 then return true end
  end
  return false
end

local addPrime = function (self, c)
  self.prime_count = self.prime_count + 1
  self.primes[self.prime_count] = c
end

local p = { prime_count=0, primes={} }
local c = 1
while p.prime_count < 1e6 do
  if not isPrimeDivisible(p, c) then
    addPrime(p, c)
  end
  c = c + 1
end
print(p.primes[p.prime_count])

--[[
Javascript code taken from V8 benchmark
V8 --prof : 9.355s
LuaJIT    : 4.604s (tables only)
MAD object: 5.499s (functions or methods)
MAD object: 5.263s (immediate methods)

function Primes() {
  this.prime_count = 0;
  this.primes = new Array(1000000);
  this.getPrimeCount = function() { return this.prime_count; }
  this.getPrime = function(i) { return this.primes[i]; }
  this.addPrime = function(i) {
    this.primes[this.prime_count++] = i;
  }

  this.isPrimeDivisible = function(candidate) {
    for (var i = 1; i < this.prime_count; ++i) {
      var current_prime = this.primes[i];
      if (current_prime * current_prime > candidate) {
        return false;
      }
      if ((candidate % this.primes[i]) == 0) return true;
    }
    return false;
  }
};

function main() {
  p = new Primes();
  var c = 1;
  while (p.getPrimeCount() < 1000000) {
    if (!p.isPrimeDivisible(c)) {
      p.addPrime(c);
    }
    c++;
  }
  print(p.getPrime(p.getPrimeCount()-1));
}

main();
]]