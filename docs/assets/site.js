(function () {
  "use strict";

  /* ---------------- Theme: light / dark / auto ---------------- */
  var root = document.documentElement;
  var STORAGE_KEY = "charm-theme";

  function applyTheme(mode) {
    if (mode === "light" || mode === "dark") {
      root.setAttribute("data-theme", mode);
    } else {
      root.removeAttribute("data-theme"); // "auto" -> follow prefers-color-scheme
    }
    document.querySelectorAll(".theme-toggle button").forEach(function (btn) {
      btn.classList.toggle("active", btn.dataset.mode === mode);
    });
  }

  function getStoredTheme() {
    try {
      return localStorage.getItem(STORAGE_KEY) || "auto";
    } catch (e) {
      return "auto";
    }
  }

  function setStoredTheme(mode) {
    try {
      localStorage.setItem(STORAGE_KEY, mode);
    } catch (e) {
      /* ignore - private browsing etc. */
    }
  }

  applyTheme(getStoredTheme());

  document.addEventListener("click", function (e) {
    var btn = e.target.closest(".theme-toggle button");
    if (!btn) return;
    var mode = btn.dataset.mode;
    setStoredTheme(mode);
    applyTheme(mode);
  });

  /* ---------------- Mobile nav menu ---------------- */
  var navToggle = document.getElementById("navmenu-toggle");
  var navPanel = document.getElementById("navmenu-panel");
  if (navToggle && navPanel) {
    navToggle.addEventListener("click", function () {
      navPanel.classList.toggle("open");
    });
    navPanel.querySelectorAll("a").forEach(function (a) {
      a.addEventListener("click", function () {
        navPanel.classList.remove("open");
      });
    });
  }

  /* ---------------- "On this page" scrollspy ---------------- */
  var pagetoc = document.getElementById("pagetoc");
  if (pagetoc) {
    var links = Array.prototype.slice.call(pagetoc.querySelectorAll("a[href^='#']"));
    var targets = links
      .map(function (a) {
        return document.querySelector(a.getAttribute("href"));
      })
      .filter(Boolean);

    function onScroll() {
      var y = window.scrollY + 100;
      var current = null;
      targets.forEach(function (t) {
        if (t.offsetTop <= y) current = t;
      });
      links.forEach(function (a) {
        a.classList.remove("active");
      });
      if (current) {
        var link = pagetoc.querySelector('a[href="#' + current.id + '"]');
        if (link) link.classList.add("active");
      }
    }
    document.addEventListener("scroll", onScroll, { passive: true });
    onScroll();
  }
})();
