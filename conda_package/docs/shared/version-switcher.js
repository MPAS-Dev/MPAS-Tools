(async function () {
    const container = document.getElementById("version-switcher");
    if (!container) return;

    const metaVersion = document.querySelector('meta[name="doc-version"]');
    const currentVersion = metaVersion ? metaVersion.content : "unknown";
    console.log("Stored current version:", currentVersion);

    async function fetchVersions() {
      try {
        const res = await fetch("/shared/versions.json");
        if (!res.ok) throw new Error();
        return await res.json();
      } catch {
        const res = await fetch("../shared/versions.json");
        return await res.json();
      }
    }

    const versions = await fetchVersions();

    const select = document.createElement("select");
    select.style.marginLeft = "1em";
    select.onchange = () => {
      window.location.href = select.value;
    };

    versions.forEach(({ version, url }) => {
      const option = document.createElement("option");
      option.value = url;
      option.textContent = version;
      if (url.includes(`/${currentVersion}/`)) {
        option.selected = true;
      }
      select.appendChild(option);
    });

    const label = document.createElement("label");
    label.textContent = "Version: ";
    label.appendChild(select);
    container.appendChild(label);
  })();
