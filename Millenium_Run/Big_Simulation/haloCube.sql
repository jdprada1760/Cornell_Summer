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
FROM    (
        SELECT *, ROW_NUMBER() OVER ( PARTITION BY f.haloID ORDER BY f.SIGMA ASC ) AS RN
        FROM  (
              --Gets the haloes and galaxies IDs that are near (nb boxes close)
              SELECT (
                         POWER(
                               (
                                  (D.x - @x)*(M.x - @x)
                                + (D.y - @y)*(M.y - @y)
                                + (D.z - @z)*(M.z - @z)

                              )*(M.x - @x)/( POWER(M.x - @x,2) + POWER(M.y - @y,2) + POWER(M.z - @z,2) )
                                - (D.x - @x)
                                , 2 )
                                +
                         POWER(
                               (
                                  (D.x - @x)*(M.x - @x)
                                + (D.y - @y)*(M.y - @y)
                                + (D.z - @z)*(M.z - @z)

                              )*(M.y - @y)/( POWER(M.x - @x,2) + POWER(M.y - @y,2) + POWER(M.z - @z,2) )
                                - (D.y - @y)
                                , 2 )
                                +
                         POWER(
                               (
                                  (D.x - @x)*(M.x - @x)
                                + (D.y - @y)*(M.y - @y)
                                + (D.z - @z)*(M.z - @z)

                              )*(M.z - @z)/( POWER(M.x - @x,2) + POWER(M.y - @y,2) + POWER(M.z - @z,2) )
                                - (D.z - @z)
                                , 2 )

                         ) AS SIGMA,
                         M.np AS NP,
                         M.haloID AS haloID
              FROM
                   (
                    SELECT Ma.haloID, Ma.x, Ma.y, Ma.z, Ma.velX, Ma.velY, Ma.velZ, Ma.np
                    FROM MPAHaloTrees..MHalo Ma
                    WHERE Ma.snapNum = @snapnum
                          AND Ma.x >= @limxh AND Ma.y >= @limyh AND Ma.z >= @limzh
                          AND Ma.x < @limMxh AND Ma.y < @limMyh AND Ma.z < @limMzh
                     ) M,
                   (
                     SELECT a.*
                     FROM
                           (SELECT D.galaxyID, D.ix, D.iy, D.iz, D.x, D.y, D.z, D.velX, D.velY, D.velZ
                           FROM MPAGalaxies..DeLucia2006a D
                           WHERE  D.snapNum = @snapnum
                                 AND D.ix >= FLOOR(@limxh/10) - @nb AND D.iy >= FLOOR(@limyh/10) - @nb AND D.iz >= FLOOR(@limzh/10) - @nb
                                 AND D.ix <= FLOOR(@limMxh/10) + @nb AND D.iy <= FLOOR(@limMyh/10) + @nb AND D.iz <= FLOOR(@limMzh/10) + @nb
                           ) a
                     INNER JOIN (SELECT L.galaxyID
                                FROM MPAGalaxies..Delucia2006a_SDSS2MASS L
                                WHERE L.snapNum = @snapnum AND L.r_sdss < -19
                              ) S
                     ON S.galaxyID = a.galaxyID
                   ) D
              WHERE --Selects only galaxies that fulfill the velocity cut
                    ABS(
                         ((@H*(D.x - M.x)/@lh + D.velX - M.velX)*(M.x - @x)+
                         (@H*(D.y - M.y)/@lh + D.velY - M.velY)*(M.y - @y) +
                         (@H*(D.z - M.z)/@lh + D.velZ - M.velZ)*(M.z - @z))/SQRT( POWER(M.x - @x,2) + POWER(M.y - @y,2) + POWER(M.z - @z,2) )

                           )  < 500
            ) f
      ) g
WHERE RN = @neigh
