function CodeBlock(el)
	local include_file = string.match(el.text, "^@include%s+(.+)$")
	if include_file then
		local file = io.open(include_file, "r")
		local content = file:read("*all")
		file:close()
		return pandoc.CodeBlock(content, el.attr)
	end
end

