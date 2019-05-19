<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">
<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">
<title>gtscript documentation</title>
<link rel="stylesheet" type="text/css" href="style.css">
</head>
<body>
<div id="menu">
<ul>
<li><a href="index.html">Overview</a></li>
<li><a href="pub/">Download</a></li>
<li><a href="https://github.com/genometools/genometools">Browse source</a></li>
<li><a href="http://github.com/genometools/genometools/issues/">Issue tracker</a></li>
<li><a href="documentation.html">Documentation</a></li>
  <ul class="submenu">
    <li><a href="tools.html">Tools</a></li>
    <li><a href="manuals.html">Manuals</a></li>
    <li><a href="libgenometools.html">C API</a></li>
    <li><a id="current" href="docs.html"><tt>gtscript</tt> docs</a></li>
    <li><a href="contract.html">Development Contract</a></li>
    <li><a href="contribute.html">Contribute</a></li>
  </ul>
<li><a href="annotationsketch.html"><tt>AnnotationSketch</tt></a></li>
<li><a href="/cgi-bin/gff3validator.cgi">GFF3 validator</a></li>
<li><a href="license.html">License</a></li>
</ul>
</div>
<div id="main">
<h1><tt>gtscript</tt> documentation</h1>
<p>
This document describes the <tt>gtscript</tt> API. <tt>gtscript</tt> is
basically <a href="http://www.lua.org/">Lua</a> (an embeddable scripting
language) plus parts of the <i>GenomeTools</i> C libraries exported to
Lua. Because the <i>GenomeTools</i> binary <tt>gt</tt> contains an embedded Lua
interpreter, <tt>gtscript</tt> files can be executed with the <tt>gt</tt>
binary.
The parts of the <i>GenomeTools</i> C libraries exported to Lua (the
<tt>gtscript</tt> API) are described here, for a documentation of Lua itself and
its APIs, please refer to the <a href="http://www.lua.org/manual/5.1/">Lua
reference manual</a>.
</p>
<p>
See the <a href="#Index">index</a> for an alphabetical list of all available
interfaces.
</p>
<h2>Notes</h2>
<ul>
<li>
You have to add <code>require 'gtlua'</code> to your script
in order to load the parts of the <tt>gtscript</tt> API which are implemented in
<tt>gtscript</tt> itself.
<li>
By default, all functions of the <tt>gtscript</tt> API are contained in the
<code>gt</code> table.  That is, you have to prepend all functions calls with
``<code>gt.</code>''. For example, write <code>gt.show()</code> to call the
<code>show()</code> function documented below.  You can use
<code>gt.export()</code> or <code>gt.re()</code> to export the <code>gt</code>
table to the global environment, which makes the prepending unnecessary (but
clutters your global environment).
<li>
If a function is documented as ``returns an array'', this means that the function
returns a table where only consecutive
integer keys from 1 on are used and which can be traversed with the
<a
href="http://www.lua.org/manual/5.1/manual.html#pdf-ipairs"><code>ipairs()</code></a> function.
</ul>
