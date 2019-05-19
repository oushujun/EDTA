<a name="<%= prototype%>"></a>
<?lua
  if args then
    func = name .. "(" .. args ..")"
  else
    func = name
  end
?>
<code><%=func%></code>
<p>
<%=comment%>
</p>
<hr>
