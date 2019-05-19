--
-- Lua command line option parser.
-- Interface based on Pythons optparse.
-- http://docs.python.org/lib/module-optparse.html
-- (c) 2008 David Manura, Licensed under the same terms as Lua (MIT license)
--
-- To be used like this:
-- t={usage="<some usage message>", version="<version string>"}
-- op = xlua.OptionParser(t)
-- op:option{"<opt>", action=<action>, dest=<dest>, help="<help message for this option>"}
--
-- with :
--   <opt> the option string to be used (can be anything, if one letter opt, then should be -x val, more letters: -xy=val )
--   <action> one of
--   - store: store in options as key, val
--   - store_true: stores key, true
--   - store_false: stores key, false
--   <dest> is the key under which the option is saved
--
-- options,args = op.parse_args()
--
-- now options is the table of options (key, val) and args is the table with non-option arguments.
-- You can use op.fail(message) for failing and op.help() for printing the usage as you like.
--
-- modifed by Benoit Corda, Clement Farabet
--

OptionParser = {}

function OptionParser:new(t)
   local self = {}
   self.usage = t.usage
   self.oneliner = t.oneliner
   self.option_descriptions = {}
   self.option_of = {}
   for k,v in pairs(OptionParser) do
      self[k] = v
   end
   self:option{"-h", "--help", action="store_true", dest="help",
               help="show this help message and exit"}
   return self
end

function OptionParser:fail(s) -- extension
   io.stderr:write(s .. '\n')
   self:help()
   os.exit(1)
end

function OptionParser:option(optdesc)
   self.option_descriptions[#self.option_descriptions+1] = optdesc
   for _,v in ipairs(optdesc) do
      self.option_of[v] = optdesc
   end
end

function OptionParser:parse(options)
   local options = options or {}
   local args = {}

   -- set the default
   for _,v in ipairs(self.option_descriptions) do
      if v.default ~= nil and options[v.dest]==nil then
         options[v.dest] = v.default
      end
   end

   if not arg then
      options.__main__ = false -- python like main
      self.options = options
      return options, args
   end
   options.__main__ = true -- python like main

   -- expand options (e.g. "--input=file" -> "--input", "file")
   local arg = {unpack(arg)}
   for i=#arg,1,-1 do local v = arg[i]
      local flag, val = v:match('^(%-%-%w+)=(.*)')
      if flag then
         arg[i] = flag
         table.insert(arg, i+1, val)
      end
   end

   local i = 1
   while i <= #arg do
      local v = arg[i]
      local optdesc = self.option_of[v]
      if optdesc then
         local default = optdesc.default
         local action = optdesc.action
         local val = default
         if action == 'store' or action == nil then
            i = i + 1
            val = arg[i] or default
            if not val then self:fail('option requires an argument ' .. v) end
         elseif action == 'store_true' then
            val = true
         elseif action == 'store_false' then
            val = false
         end
         options[optdesc.dest] = val
      else
         if v:match('^%-') then self:fail('invalid option ' .. v) end
         args[#args+1] = v
      end
      i = i + 1
   end
   for k,opt in pairs(self.option_of) do
      if opt.req and not options[opt.dest] then
         self:fail('option '.. k .. ' requires an argument ')
      end
   end
   if options.help then
      self:help()
      os.exit()
   end
   -- set the default if nil
   self.options = options
   return options, args
end

function OptionParser:flags(optdesc)
   local sflags = {}
   local action = optdesc and optdesc.action
   for _,flag in ipairs(optdesc) do
      local sflagend
      if action == nil or action == 'store' then
         local metavar = optdesc.metavar or optdesc.dest:upper()
         sflagend = #flag == 2 and ' ' .. metavar
            or  '=' .. metavar
      else
         sflagend = ''
      end
      sflags[#sflags+1] = flag .. sflagend
   end
   return table.concat(sflags, ', ')
end

function OptionParser:help()
   io.stdout:write(self.oneliner .. "\n")
   if arg[-1] then
      io.stdout:write("Usage: " .. self.usage:gsub('%%prog', (arg[-1] .. ' ' .. arg[0])) .. "\n")
   elseif arg[0] then
      io.stdout:write("Usage: " .. self.usage:gsub('%%prog', arg[0]) .. "\n")
   else
      io.stdout:write("Usage: " .. self.usage:gsub('%%prog', 'THISPROG') .. "\n")
   end
   io.stdout:write("\n")
   io.stdout:write("Options:\n")
   pad = 0
   for _,optdesc in ipairs(self.option_descriptions) do
      pad = math.max(pad, #self:flags(optdesc))
   end
   for _,optdesc in ipairs(self.option_descriptions) do
      local defstr = ''
      if optdesc.req then
         defstr = ' [REQUIRED]'
      elseif optdesc.default then
         defstr = ' [default = ' .. tostring(optdesc.default) .. ']'
      end
      io.stdout:write("  " .. self:flags(optdesc) ..
                   string.rep(' ', pad - #self:flags(optdesc)) ..
                "  " .. optdesc.help .. defstr .. "\n")
   end
end

function OptionParser:tostring(generatefilename, params)
   local str = ''
   if not generatefilename then
      str = '<'.. ((arg and arg[0]) or 'interpreted.lua'):gsub('.lua','') .. "> configuration:\n"
      for k,v in pairs(self.options) do
         str = str .. ' + ' .. k .. ' = ' .. tostring(v) .. '\n'
      end
   else
      local first = true
      for i,entry in ipairs(self.option_descriptions) do
         local key = entry[1]
         local match = true
         if #params > 0 then
            match = false
            for i,param in ipairs(params) do
               if key == param then match = true; break end
            end
         end
         local val = self.options[entry.dest]
         if val and match then
            if first then
               str = str .. key .. '=' .. tostring(val)
            else
               str = str .. ',' .. key .. '=' .. tostring(val)
            end
            first = false
         end
      end
      str = str:gsub('/','|'):gsub(' ','_')
   end
   return str
end

function OptionParser:summarize(compact)
   io.write(self:tostring(compact))
end
