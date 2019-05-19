<a name="Index"></a>
<h2>Index</h2>
<%for _, name in ipairs(names) do%>
  <a href="#<%= name%>"><code><%= name%></code></a><br>
<%end%>
