# Motor CD Shunt — control de velocidad por campo (Scilab / Xcos)

| Archivo | Contenido |
|---|---|
| [`GUIA_DE_ESTUDIO.md`](GUIA_DE_ESTUDIO.md) | **Empieza aquí.** Explicación pedagógica completa, receta para exámenes y ejercicios |
| `calculos_motor_shunt.sce` | Parte A + Parte B (ode, simulación del .zcos y las 2 gráficas). Crea `Diaz_Diego.zcos` si no existe |
| `generar_xcos.sce` | Construye el diagrama de bloques por código |
| `Diaz_Diego.zcos` | Diagrama Xcos entregable. **No viene descargado:** lo crea tu Scilab al ejecutar `calculos_motor_shunt.sce` (así es compatible con tu versión) |
| `informe/Diaz_Diego_Informe.docx` | Informe final en Word (pasar a PDF para entregar) |
| `figuras/` | Gráficas y captura del diagrama |
| `resultados_consola.txt` | Salida completa del script de Scilab |
| `verificacion_python.py` | Verificación independiente con SciPy |

Resultados clave: Ke = Kt = 0.65950 V·s/rad, Ia = 19 A, Ea = 110.5 V, **Rf' = 136.63 Ω**.

```scilab
exec("calculos_motor_shunt.sce", -1)
```
Probado en Scilab 2024.0.0.
