// MathJax (version 2, as loaded by FORD) settings for the API docs: inline math
// between single dollar signs too, as GitHub renders it (the README is included
// in the docs' front page). FORD loads this file before MathJax itself, so it
// sets the configuration object that MathJax reads when it starts.
window.MathJax = {
  tex2jax: {
    inlineMath: [['$', '$'], ['\\(', '\\)']],
    displayMath: [['$$', '$$'], ['\\[', '\\]']],
    processEscapes: true
  }
};
