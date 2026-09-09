// ============================================
// FACTOR R - EROSIVIDAD DE LA LLUVIA
// Fuente: CHIRPS Daily (~5.5 km resolución nativa)
// Metodología: Índice Modificado de Fournier (IMF) -> R (Arnoldus, 1980)
// Período: 2000-2025 (26 años, promediados)
// Departamento: SACATEPÉQUEZ
// ============================================

var roi = ee.FeatureCollection('projects/fondecyt-entorno-urbano/assets/sacatepequez_limite');

var chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY').filterBounds(roi);

var calcularRAnual = function(year) {
  year = ee.Number(year);
  var meses = ee.List.sequence(1, 12);

  var precipMensual = ee.ImageCollection.fromImages(
    meses.map(function(m) {
      var inicio = ee.Date.fromYMD(year, m, 1);
      var fin = inicio.advance(1, 'month');
      return chirps.filterDate(inicio, fin).sum().rename('precip_mes').set('mes', m);
    })
  );

  var precipAnual = precipMensual.sum().rename('precip_anual');

  // IMF = suma( (precip_mes^2) / precip_anual )
  var imfMensual = precipMensual.map(function(img) {
    return img.pow(2).divide(precipAnual).rename('imf_mes');
  });
  var imfAnual = imfMensual.sum().rename('IMF');

  // R = 2.56 * IMF^1.065 (Arnoldus, 1980)
  return imfAnual.pow(1.065).multiply(2.56).rename('R').set('year', year);
};

var anios = ee.List.sequence(2000, 2025);
var rPorAnio = ee.ImageCollection.fromImages(anios.map(calcularRAnual));
var factorR = rPorAnio.mean().rename('R_factor').clip(roi);

// ============================================
// VALIDACIÓN — R_mean: 1015.2 | rango: 609.2 - 1737.6
// ============================================
Map.centerObject(roi, 10);
Map.addLayer(factorR, {min: 0, max: 800, palette: ['blue','green','yellow','red']}, 'Factor R');

print('Estadísticas Factor R (Sacatepéquez):', factorR.reduceRegion({
  reducer: ee.Reducer.minMax().combine(ee.Reducer.mean(), '', true),
  geometry: roi.geometry(),
  scale: 5566,
  maxPixels: 1e9
}));

// ============================================
// FACTOR LS - LONGITUD Y GRADO DE PENDIENTE
// DEM: Copernicus DEM GLO-30 (30m) | Flow Accumulation: HydroSHEDS 15ACC (~464m)
// Componente L: enfoque de acumulación de flujo (Moore & Burch, 1986), exponente
//   fijo m=0.4 -- simplificación práctica; el valor rigurosamente varía 0.2-0.5
//   según pendiente (Foster & Wischmeier)
// Componente S: McCool et al. (1987)
// Límite de 122m en longitud de flujo, siguiendo recomendación de Fu et al. (2006),
//   Yang (2015) y Borrelli et al. (2017) -- no atribuible a Desmet & Govers (1996)
// Flow accumulation remuestreado de ~464m a 30m vía interpolación bilineal para
//   compatibilizar con el DEM -- simplificación sobre variable acumulativa, sin
//   cuantificar su efecto en la precisión
// Departamento: SACATEPÉQUEZ
// ============================================
// El .mosaic() no hereda la proyección nativa; sin reasignarla, ee.Terrain.slope() falla
var demColeccion = ee.ImageCollection('COPERNICUS/DEM/GLO30').filterBounds(roi).select('DEM');
var proyeccionNativa = demColeccion.first().projection();
var dem = demColeccion.mosaic().setDefaultProjection(proyeccionNativa).clip(roi);

var elevacion = dem.rename('elevacion_m');
var pendienteGrados = ee.Terrain.slope(dem);

var flowAcc = ee.Image('WWF/HydroSHEDS/15ACC').select('b1');
var flowAccResampled = flowAcc.resample('bilinear').reproject({crs: proyeccionNativa, scale: 30}).clip(roi);

var tamañoCeldaFlowAcc = 463.83; // resolución nativa HydroSHEDS 15ACC

var longitudFlujo = flowAccResampled.multiply(tamañoCeldaFlowAcc);
var longitudFlujoLimitada = longitudFlujo.min(122);
var factorL = longitudFlujoLimitada.divide(22.13).pow(0.4);

var factorS = pendienteGrados.expression(
  "slope < 9 ? 10.8 * sin(slope * pi / 180) + 0.03 : 16.8 * sin(slope * pi / 180) - 0.5",
  {'slope': pendienteGrados, 'pi': Math.PI}
);

var factorLS = factorL.multiply(factorS).rename('LS_factor');

// ============================================
// VALIDACIÓN — Elevación: 644-3970m | Pendiente: 0-65.6° | LS: 0.06-29.3 (mean 8.5)
// ============================================
Map.addLayer(elevacion, {min: 1500, max: 3500, palette: ['blue','green','yellow','brown','white']}, 'Elevación');
Map.addLayer(pendienteGrados, {min: 0, max: 65, palette: ['green','yellow','red']}, 'Pendiente (grados)');
Map.addLayer(factorLS, {min: 0, max: 30, palette: ['green','yellow','red']}, 'Factor LS');

print('Estadísticas Elevación:', elevacion.reduceRegion({
  reducer: ee.Reducer.minMax().combine(ee.Reducer.mean(), '', true),
  geometry: roi.geometry(), scale: 30, maxPixels: 1e9
}));
print('Estadísticas Pendiente:', pendienteGrados.reduceRegion({
  reducer: ee.Reducer.minMax().combine(ee.Reducer.mean(), '', true),
  geometry: roi.geometry(), scale: 30, maxPixels: 1e9
}));
print('Estadísticas Factor LS:', factorLS.reduceRegion({
  reducer: ee.Reducer.minMax().combine(ee.Reducer.mean(), '', true),
  geometry: roi.geometry(), scale: 30, maxPixels: 1e9
}));

// ============================================
// FACTOR K - ERODABILIDAD DEL SUELO
// Fuente: SoilGrids v2.0 (ISRIC), 250m resolución
// Metodología: Williams (1995), vía SWAT+ theoretical documentation
// K_USLE = f_csand * f_cl-si * f_orgc * f_hisand, convertido a SI (x0.1317, Foster 1981)
// https://swatplus.gitbook.io/io-docs/theoretical-documentation/section-4-erosion/sediment/musle/soil-erodibility-factor
// Departamento: SACATEPÉQUEZ
// ============================================

var arena   = ee.Image('projects/soilgrids-isric/sand_mean').select('sand_0-5cm_mean').divide(10).clip(roi);
var arcilla = ee.Image('projects/soilgrids-isric/clay_mean').select('clay_0-5cm_mean').divide(10).clip(roi);
var limo    = ee.Image('projects/soilgrids-isric/silt_mean').select('silt_0-5cm_mean').divide(10).clip(roi);
var soc_raw = ee.Image('projects/soilgrids-isric/soc_mean').select('soc_0-5cm_mean').clip(roi); // dg/kg

var orgC = soc_raw.divide(100); // dg/kg -> %

var fcsand = arena.multiply(-0.256).multiply(
    ee.Image(1).subtract(limo.divide(100))
  ).exp().multiply(0.3).add(0.2);

var fclsi = limo.divide(arcilla.add(limo)).pow(0.3);

var forgc = ee.Image(1).subtract(
  orgC.multiply(0.25).divide(
    orgC.add(ee.Image(3.72).subtract(orgC.multiply(2.95)).exp())
  )
);

var unoMenosArena = ee.Image(1).subtract(arena.divide(100));
var fhisand = ee.Image(1).subtract(
  unoMenosArena.multiply(0.7).divide(
    unoMenosArena.add(ee.Image(-5.51).add(unoMenosArena.multiply(22.9)).exp())
  )
);

var factorK_US = fcsand.multiply(fclsi).multiply(forgc).multiply(fhisand);
var factorK = factorK_US.multiply(0.1317).rename('K_factor');

// Relleno de huecos: SoilGrids presentó cobertura incompleta en ~8.5% del área
// (8,146 de 8,905 píxeles válidos a 250m). Se rellena con el promedio del área
// válida -- imputación simple, documentada como limitación metodológica.
var promedioK = factorK.reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: roi.geometry(), scale: 250, maxPixels: 1e9, bestEffort: true
});
var factorK_relleno = factorK.unmask(ee.Number(promedioK.get('K_factor'))).rename('K_factor');

// ============================================
// VALIDACIÓN — K: 0.0147-0.0178 (mean 0.0165) | Cobertura SoilGrids: 91.5%
// ============================================
Map.addLayer(factorK_relleno, {min: 0, max: 0.07, palette: ['green','yellow','red']}, 'Factor K (relleno)');

print('Promedio K usado para relleno:', promedioK);
print('Estadísticas Factor K (con relleno):', factorK_relleno.reduceRegion({
  reducer: ee.Reducer.minMax().combine(ee.Reducer.mean(), '', true),
  geometry: roi.geometry(), scale: 250, maxPixels: 1e9, bestEffort: true
}));
print('Píxeles válidos K (departamento completo):', factorK.reduceRegion({
  reducer: ee.Reducer.count(),
  geometry: roi.geometry(), scale: 250, maxPixels: 1e9, bestEffort: true
}));
print('Píxeles totales esperados (departamento completo):', ee.Image(1).reduceRegion({
  reducer: ee.Reducer.count(),
  geometry: roi.geometry(), scale: 250, maxPixels: 1e9, bestEffort: true
}));

print('Percentiles Factor K:', factorK_relleno.reduceRegion({
  reducer: ee.Reducer.percentile([5, 25, 50, 75, 95]),
  geometry: roi.geometry(),
  scale: 250,
  maxPixels: 1e9,
  bestEffort: true
}));

// ============================================
// FACTOR C - COBERTURA VEGETAL
// Fuente: Sentinel-2 SR (10m), período 2015-2025
// (Limitación: Sentinel-2 no cubre 2000-2014; ver introducción de la bitácora)
// Metodología: NDVI -> Factor C, Van der Knijff, Jones & Montanarella (2000)
//   C = exp[-α·NDVI/(β-NDVI)], con α=2, β=1 (parámetros estándar verificados)
// Limitación conocida de la fórmula: puede sobreestimar C (subestimar protección
//   vegetal) cuando NDVI < 0.65 (Fiorucci et al., D-RUSLE)
// C se calculó a partir del NDVI promediado del período 2015-2025, no promediando
//   C por imagen individual -- simplificación práctica no verificada contra el
//   procedimiento original de Van der Knijff et al.
// Correcciones aplicadas:
//   1. NDVI clampeado a mínimo 0.01 -- evita que NDVI negativo (agua/sombras)
//      rompa el rango teórico de C [0,1] vía exponente positivo
//   2. Agua excluida vía JRC Global Surface Water (occurrence > 50%), con
//      unmask(0) por enmascarado nativo del dataset en tierra firme
// Departamento: SACATEPÉQUEZ
// ============================================

var sentinel2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(roi)
  .filterDate('2015-06-01', '2025-12-31')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20));

function enmascararNubes(imagen) {
  var scl = imagen.select('SCL');
  var mascara = scl.neq(3).and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10));
  return imagen.updateMask(mascara);
}
var sentinel2Limpio = sentinel2.map(enmascararNubes);

function calcularNDVI(imagen) {
  return imagen.addBands(imagen.normalizedDifference(['B8', 'B4']).rename('NDVI'));
}
var conNDVI = sentinel2Limpio.map(calcularNDVI);

var ndviPromedio = conNDVI.select('NDVI').mean().clip(roi);
var ndviClampeado = ndviPromedio.max(0.01);

var aguaRaw = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('occurrence');
var esAgua = aguaRaw.unmask(0).gt(50);

var factorC = ndviClampeado.expression("exp(-2 * ndvi / (1 - ndvi))", {'ndvi': ndviClampeado})
  .min(1).max(0)
  .updateMask(esAgua.not())
  .rename('C_factor');

// ============================================
// VALIDACIÓN (departamento completo)
// NDVI  -- mean: 0.613 | P5: 0.252 | P50: 0.670 | P95: 0.807 | max: 0.882
// C     -- mean: 0.105 | P5: 0.0013 | P50: 0.0175 | P95: 0.514 | max: 0.980
// Cobertura: 5,572,377 píxeles válidos (agua excluida correctamente)
// ============================================

Map.addLayer(ndviPromedio, {min: -0.2, max: 0.9, palette: ['red','yellow','green']}, 'NDVI promedio');
Map.addLayer(factorC, {min: 0, max: 1, palette: ['green','yellow','red']}, 'Factor C');

print('Percentiles NDVI (Sacatepéquez completo):', ndviPromedio.reduceRegion({
  reducer: ee.Reducer.percentile([5, 25, 50, 75, 95]).combine(ee.Reducer.minMax(), '', true).combine(ee.Reducer.mean(), '', true),
  geometry: roi.geometry(), scale: 10, maxPixels: 1e9, bestEffort: true, tileScale: 4
}));

print('Percentiles Factor C (Sacatepéquez completo):', factorC.reduceRegion({
  reducer: ee.Reducer.percentile([5, 25, 50, 75, 95]).combine(ee.Reducer.minMax(), '', true).combine(ee.Reducer.mean(), '', true),
  geometry: roi.geometry(), scale: 10, maxPixels: 1e9, bestEffort: true, tileScale: 4
}));

print('Píxeles válidos Factor C (departamento completo):', factorC.reduceRegion({
  reducer: ee.Reducer.count(),
  geometry: roi.geometry(), scale: 10, maxPixels: 1e9, bestEffort: true, tileScale: 4
}));

// ============================================
// FACTOR P - PRÁCTICAS DE CONSERVACIÓN
// Sin datos de campo disponibles sobre prácticas de conservación (terrazas,
// curvas de nivel, etc.) -> se asume P = 1 (sin prácticas), decisión
// metodológica acordada con el cliente
// Departamento: SACATEPÉQUEZ
// ============================================

var factorP = ee.Image(1).clip(roi).rename('P_factor');

// ============================================
// RESULTADO FINAL - ECUACIÓN RUSLE
// A = R x K x LS x C x P
// ============================================

var erosionRUSLE = factorR
  .multiply(factorK_relleno)
  .multiply(factorLS)
  .multiply(factorC)
  .multiply(factorP)
  .rename('erosion_ton_ha_year');

Map.addLayer(erosionRUSLE, {min: 0, max: 200, palette: ['green','yellow','orange','red','darkred']}, 'Erosión RUSLE (ton/ha/año)');

// ============================================
// VALIDACIÓN (departamento completo)
// ============================================

print('Percentiles Erosión RUSLE (Sacatepéquez completo):', erosionRUSLE.reduceRegion({
  reducer: ee.Reducer.percentile([5, 25, 50, 75, 95]),
  geometry: roi.geometry(), 
  scale: 30,  // <- bajamos de 10 a 30 solo para esta validación estadística
  maxPixels: 1e9, 
  bestEffort: true, 
  tileScale: 8  // <- subimos el tileScale para repartir más el cómputo
}));

// ============================================
// EXPORTACIÓN A GOOGLE DRIVE
// ============================================

Export.image.toDrive({
  image: erosionRUSLE,
  description: 'RUSLE_erosion_sacatepequez_final',
  folder: 'Cartos_Erosion_GEE',
  region: roi.geometry(), scale: 10, crs: 'EPSG:4326', maxPixels: 1e10
});
Export.image.toDrive({
  image: factorR.rename('R'),
  description: 'RUSLE_factorR_sacatepequez',
  folder: 'Cartos_Erosion_GEE',
  region: roi.geometry(), scale: 10, crs: 'EPSG:4326', maxPixels: 1e10
});
Export.image.toDrive({
  image: factorK_relleno.rename('K'),
  description: 'RUSLE_factorK_sacatepequez_final',
  folder: 'Cartos_Erosion_GEE',
  region: roi.geometry(), scale: 10, crs: 'EPSG:4326', maxPixels: 1e10
});
Export.image.toDrive({
  image: factorLS.rename('LS'),
  description: 'RUSLE_factorLS_sacatepequez',
  folder: 'Cartos_Erosion_GEE',
  region: roi.geometry(), scale: 10, crs: 'EPSG:4326', maxPixels: 1e10
});
Export.image.toDrive({
  image: factorC.rename('C'),
  description: 'RUSLE_factorC_sacatepequez',
  folder: 'Cartos_Erosion_GEE',
  region: roi.geometry(), scale: 10, crs: 'EPSG:4326', maxPixels: 1e10
});
Export.image.toDrive({
  image: elevacion,
  description: 'DEM_elevacion_sacatepequez',
  folder: 'Cartos_Erosion_GEE',
  region: roi.geometry(), scale: 10, crs: 'EPSG:4326', maxPixels: 1e10
});

