--Position of the point of view
DECLARE @x FLOAT
DECLARE @y FLOAT
DECLARE @z FLOAT
--Snapshot to be analized
DECLARE @snapnum INT
--Hubble constant
DECLARE @H FLOAT
--Neighbour
DECLARE @neigh INT

--Asignment of variables
SET @snapnum = 63
SET @neigh = 4
SET @x = 250
SET @y = 250
SET @z = 250
SET @H = 100

            (SELECT SQRT
                   POWER(  ((c.gx - @x)*c.dx + (c.gy - @y)*c.dy + (c.gz - @z)*c.dz)*c.dx - (c.gx - @x)  , 2 )  +
                   POWER(  ((c.gx - @x)*c.dx + (c.gy - @y)*c.dy + (c.gz - @z)*c.dz)*c.dy - (c.gy - @y)  , 2 )  +
                   POWER(  ((c.gx - @x)*c.dx + (c.gy - @y)*c.dy + (c.gz - @z)*c.dz)*c.dz - (c.gz - @z)  , 2 )
                              ) AS DIST,
             FROM
                  -- Selects the position, velocities and position, velocities of galaxies and haloes and the unit vectors to haloes
                  (SELECT D.x AS gx,
                          D.y AS gy,
                          D.z AS gz,
                          M.x AS hx,
                          M.y AS hy,
                          M.z AS hz,
                          D.velX AS gvelX,
                          D.velY AS gvelY,
                          D.velZ AS gvelZ,
                          M.velX AS hvelX,
                          M.velY AS hvelY,
                          M.velZ AS hvelZ,
                          (M.x - @x)/SQRT( POWER(M.x - @x,2) + POWER(M.y - @y,2) + POWER(M.z - @z,2) ) AS dx,
                          (M.y - @y)/SQRT( POWER(M.x - @x,2) + POWER(M.y - @y,2) + POWER(M.z - @z,2) ) AS dy,
                          (M.z - @z)/SQRT( POWER(M.x - @x,2) + POWER(M.y - @y,2) + POWER(M.z - @z,2) ) AS dz,
                   FROM DeLucia2006a D,
                        MPAHalo M
                   WHERE --Selects snapshot
                         D.snapnum = @snapnum
                         AND M.snapNum = @snapnum
                         --Omits galaxies with r abs magnitude < -19
                         AND D.mag_r < -19
                  )


             WHERE
                  --Only counts galxies with relative velocities to the halo < 500 Km/s
                  ABS(
                  --Radial component of the total velocity of the galaxy
                      c.dx*(@H*(c.gx-@x) + c.gvelX) + c.dy*(@H*(c.gy-@y) + c.gvelY) + c.dz*(@H*(c.gz-@z) + c.gvelZ)
                  --Radial component of the total velocity of the halo
                      -(c.dx*(@H*(c.hx-@x) + c.hvelX) + c.dy*(@H*(c.hy-@y) + c.hvelY) + c.dz*(@H*(c.hz-@z) + c.hvelZ))

                     ) < 500
             ) a
