--Position of the point of view
DECLARE @x FLOAT
DECLARE @y FLOAT
DECLARE @z FLOAT
--Hubble constant
DECLARE @H FLOAT
-- Little h
DECLARE @lh float
--Neighbour
DECLARE @neigh INT
--Limits the box of halos (this time uses the number asociated to each position pointer)
DECLARE @limxh INT
DECLARE @limyh INT
DECLARE @limzh INT
DECLARE @limMxh INT
DECLARE @limMyh INT
DECLARE @limMzh INT
--Snapshot to be analyzed
DECLARE @snapnum INT
--Limits the size of the cube of galaxies for each halo
DECLARE @nb INT

--Asignment of variables
SET @nb = 3
SET @snapnum = 63
SET @limxh = %x
SET @limyh = %y
SET @limzh = %z
SET @limMxh = $x
SET @limMyh = $y
SET @limMzh = $z
SET @lh = 0.7
SET @snapnum = 63
SET @neigh = 4
SET @x = 0
SET @y = 0
SET @z = 0
SET @H = 100

SELECT *
FROM
      (SELECT  f.DIST,
               f.NP,
               ROW_NUMBER() OVER ( PARTITION BY f.haloID ORDER BY f.DIST ASC ) AS RN
      FROM
            (SELECT
                  POWER(  ((c.gx - @x)*c.dx + (c.gy - @y)*c.dy + (c.gz - @z)*c.dz)*c.dx - (c.gx - @x)  , 2 )  +
                  POWER(  ((c.gx - @x)*c.dx + (c.gy - @y)*c.dy + (c.gz - @z)*c.dz)*c.dy - (c.gy - @y)  , 2 )  +
                  POWER(  ((c.gx - @x)*c.dx + (c.gy - @y)*c.dy + (c.gz - @z)*c.dz)*c.dz - (c.gz - @z)  , 2 )
                  AS DIST,
                  c.haloID AS haloID,
                  c.NP AS NP
            FROM
                -- Selects the position, velocities and position, velocities of galaxies and haloes and the unit vectors to haloes
                (SELECT M.haloID AS haloID,
                        M.np AS NP,
                        D.x/@lh AS gx,
                        D.y/@lh AS gy,
                        D.z/@lh AS gz,
                        M.x/@lh AS hx,
                        M.y/@lh AS hy,
                        M.z/@lh AS hz,
                        --Difference in peculiar velocities
                        D.velX - M.velX AS dvx,
                        D.velY - M.velY AS dvy,
                        D.velZ - M.velZ AS dvz,
                        (M.x - @x)/SQRT( POWER(M.x - @x,2) + POWER(M.y - @y,2) + POWER(M.z - @z,2) ) AS dx,
                        (M.y - @y)/SQRT( POWER(M.x - @x,2) + POWER(M.y - @y,2) + POWER(M.z - @z,2) ) AS dy,
                        (M.z - @z)/SQRT( POWER(M.x - @x,2) + POWER(M.y - @y,2) + POWER(M.z - @z,2) ) AS dz
                 FROM
                       --Selects the box of galaxies (filter by snapshot)
                       (
                       SELECT a.*
                       FROM
                             (SELECT D.galaxyID, D.ix, D.iy, D.iz, D.x, D.y, D.z, D.velX, D.velY, D.velZ
                             FROM MPAGalaxies..DeLucia2006a D
                             WHERE --Selects snapshot
                                   D.snapNum = @snapnum
                                   --Restricts the box of galaxies (bigger than haloes)
                                   AND D.ix >= FLOOR(@limxh/10) - @nb AND D.iy >= FLOOR(@limyh/10) - @nb AND D.iz >= FLOOR(@limzh/10) - @nb
                                   AND D.ix <= FLOOR(@limMxh/10) + @nb AND D.iy <= FLOOR(@limMyh/10) + @nb AND D.iz <= FLOOR(@limMzh/10) + @nb
                             ) a
                       INNER JOIN (SELECT L.galaxyID
                                  FROM MPAGalaxies..Delucia2006a_SDSS2MASS L
                                  WHERE L.snapNum = @snapnum AND L.r_sdss < -19
                                ) S

                       ON S.galaxyID = a.galaxyID
                     ) D,
                      --Selects the box of haloes to work
                      (
                      SELECT Ma.haloID, Ma.x, Ma.y, Ma.z, Ma.velX, Ma.velY, Ma.velZ, Ma.np
                      FROM MPAHaloTrees..MHalo Ma
                      WHERE --Selects snapshot
                            Ma.snapNum = @snapnum
                            --Restricts the box of haloes
                            AND Ma.x >= @limxh AND Ma.y >= @limyh AND Ma.z >= @limzh
                            AND Ma.x < @limMxh AND Ma.y < @limMyh AND Ma.z < @limMzh
                       ) M
                ) c
            WHERE --Only counts galxies with relative radial velocities to the halo < 500 Km/s
                 ABS(
                     -- Radial component of the difference in peculiar velocities
                     c.dvx*c.dx + c.dvy*c.dy + c.dvz*c.dz +
                     -- Radial component of the diference in velocities due to Hubbles law
                     @H*((c.gx - c.hx)*dx + (c.gy - c.hy)*dy + (c.gz - c.hz)*dz)

                   )  < 500
           ) f
      ) g
WHERE RN = @neigh





--Selects the box of haloes to work
SELECT M.haloID, M.x, M.y, M.z, M.velX, M.velY, M.velZ, M.np
INTO HaloBox
FROM MPAHaloTrees..MHalo M
WHERE --Selects snapshot
      M.snapNum = @snapnum
      --Restricts the box of haloes
      AND M.x >= @limxh AND M.y >= @limyh AND M.z >= @limzh
      AND M.x < @limMxh AND M.y < @limMyh AND M.z < @limMzh


--Selects the box of galaxies (filter by snapshot)
SELECT a.*
FROM
      (SELECT D.galaxyID, D.ix, D.iy, D.iz, D.x, D.y, D.z, D.velX, D.velY, D.velZ
      FROM MPAGalaxies..DeLucia2006a D
      WHERE --Selects snapshot
            D.snapNum = @snapnum
            --Restricts the box of galaxies (bigger than haloes)
            AND D.ix >= FLOOR(@limxh/10) - @nb AND D.iy >= FLOOR(@limyh/10) - @nb AND D.iz >= FLOOR(@limzh/10) - @nb
            AND D.ix <= FLOOR(@limMxh/10) + @nb AND D.iy <= FLOOR(@limMyh/10) + @nb AND D.iz <= FLOOR(@limMzh/10) + @nb
      ) a
INNER JOIN (SELECT L.galaxyID
           FROM MPAGalaxies..Delucia2006a_SDSS2MASS L
           WHERE L.snapNum = @snapnum AND L.r_sdss < -19
         ) S
ON S.galaxyID = a.galaxyID
