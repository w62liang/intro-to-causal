document.addEventListener("DOMContentLoaded", function () {
  // 1) 加载 KaTeX CSS
  var l = document.createElement("link");
  l.rel = "stylesheet";
  l.href = "https://cdn.jsdelivr.net/npm/katex@0.16.9/dist/katex.min.css";
  document.head.appendChild(l);

  // 2) 加载 KaTeX 主库 + auto-render 插件
  var s1 = document.createElement("script");
  s1.src = "https://cdn.jsdelivr.net/npm/katex@0.16.9/dist/katex.min.js";
  s1.onload = function () {
    var s2 = document.createElement("script");
    s2.src = "https://cdn.jsdelivr.net/npm/katex@0.16.9/dist/contrib/auto-render.min.js";
    s2.onload = function () {
      function render() {
        renderMathInElement(document.body, {
          delimiters: [
            { left: "$$", right: "$$", display: true },
            { left: "\\[", right: "\\]", display: true },
            { left: "$",  right: "$",  display: false },
            { left: "\\(", right: "\\)", display: false }
          ],
          throwOnError: false
        });
      }
      render();

      // ⚡ 兼容 MkDocs Material 的 instant loading
      if (typeof document$ !== "undefined") {
        document$.subscribe(render);
      }
    };
    document.body.appendChild(s2);
  };
  document.body.appendChild(s1);
});

var style = document.createElement("style");
style.innerHTML = `
  .katex .mathnormal {
    font-family: "Roboto", cursive !important;
  }
`;
document.head.appendChild(style);
