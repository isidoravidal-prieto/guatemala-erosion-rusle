# Guatemala Hydric Erosion Risk Mapping (RUSLE)

Google Earth Engine script implementing the RUSLE (Revised Universal Soil Loss Equation) methodology to estimate annual soil water erosion risk across four departments in Guatemala (Alta Verapaz, Quiché, Sololá, Sacatepéquez), for an environmental engineering client.

Combines five factors natively in Earth Engine:
- **R** (rainfall erosivity) — CHIRPS daily precipitation, 2000–2025
- **K** (soil erodibility) — SoilGrids v2.0 texture and organic carbon
- **LS** (slope length and steepness) — Copernicus DEM GLO-30 + HydroSHEDS
- **C** (vegetation cover) — Sentinel-2 NDVI, 2015–2025
- **P** (conservation practices) — held at baseline, no field data available

Full write-up and results: [isidoravidal-prieto.github.io/projects/erosion-guatemala.html](https://isidoravidal-prieto.github.io/projects/erosion-guatemala.html)

---

# Mapas de Erosión Hídrica del Suelo, Guatemala (RUSLE)

Script de Google Earth Engine que implementa la metodología RUSLE (Revised Universal Soil Loss Equation) para estimar el riesgo anual de erosión hídrica del suelo en cuatro departamentos de Guatemala (Alta Verapaz, Quiché, Sololá, Sacatepéquez), para un cliente de ingeniería ambiental.

Combina cinco factores directamente en Earth Engine:
- **R** (erosividad de la lluvia) — precipitación diaria CHIRPS, 2000–2025
- **K** (erodabilidad del suelo) — textura y carbono orgánico de SoilGrids v2.0
- **LS** (longitud y grado de pendiente) — Copernicus DEM GLO-30 + HydroSHEDS
- **C** (cobertura vegetal) — NDVI de Sentinel-2, 2015–2025
- **P** (prácticas de conservación) — asumido en línea base, sin datos de campo disponibles

Informe completo y resultados: [isidoravidal-prieto.github.io/projects/erosion-guatemala.html](https://isidoravidal-prieto.github.io/projects/erosion-guatemala.html)
