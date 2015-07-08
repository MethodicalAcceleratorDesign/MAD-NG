local M = { __author = 'ldeniau', __version = '2015.06', help = {}, test = {} }

--for k,v in pairs(arg) do
--	print(k,v)
--end

-- skip luajit options for now...
while arg[1]:sub(1,1) == '-' do
	table.remove(arg,1)
end

dofile( table.remove(arg,1) )