@echo off

FOR /L %%A IN (10,1,10) DO (
	FOR /L %%B IN (10,1,13) DO (
		FOR /L %%C IN (10,1,10) DO (
			FOR /L %%D IN (10,1,14) DO (
				FOR /L %%E IN (10,1,13) DO (
					FOR /L %%F IN (10,1,11) DO (
						FOR /L %%G IN (10,1,11) DO (
							FOR /L %%H IN (10,1,10) DO (
								FOR /L %%I IN (10,1,11) DO (
									FOR /L %%J IN (10,1,12) DO (
										Scheduling.exe %%A %%B %%C %%D %%E %%F %%G %%H %%I %%J
									)
								)
							)
						)
					)
				)
			)
		)
	)
)