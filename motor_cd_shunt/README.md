# Motor CD Shunt — control de velocidad por campo (Scilab / Xcos)

| Archivo | Contenido |
|---|---|
| [`GUIA_DE_ESTUDIO.md`](GUIA_DE_ESTUDIO.md) | **Empieza aquí.** Explicación pedagógica completa, receta para exámenes y ejercicios |
| `calculos_motor_shunt.sce` | Parte A (memoria de cálculo) + Parte B (ode, simulación del .zcos, gráficas) |
| `generar_xcos.sce` | Construye el diagrama de bloques por código |
| `APELLIDO_NOMBRE.zcos` | Diagrama Xcos entregable (renómbralo) |
| `informe/` | Borrador del informe: `.md` (fuente), `.pdf`, `.docx` (editable en Word) |
| `figuras/` | Gráficas y captura del diagrama |
| `resultados_consola.txt` | Salida completa del script de Scilab |
| `verificacion_python.py` | Verificación independiente con SciPy |

Resultados clave: Ke = Kt = 0.65950 V·s/rad, Ia = 19 A, Ea = 110.5 V, **Rf' = 136.63 Ω**.

```scilab
exec("calculos_motor_shunt.sce", -1)
```
Probado en Scilab 2024.0.0.
