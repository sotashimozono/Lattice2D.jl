@testset "Index Accessor Functions" begin
    @testset "Square Lattice" begin
        Lx, Ly = 4, 3
        lat = build_lattice(Square, Lx, Ly)

        @testset "get_site_index" begin
            # Test basic indexing
            idx = get_site_index(lat, 1, 1)
            @test idx == 1

            idx = get_site_index(lat, 2, 1)
            @test idx == 2

            idx = get_site_index(lat, 1, 2)
            @test idx == 5  # (0 + (2-1)*4) * 1 + 1 = 5

            # Test all positions
            for x in 1:Lx, y in 1:Ly
                idx = get_site_index(lat, x, y)
                @test idx == lat.site_map[x, y]
            end

            # Test boundary checks
            @test_throws ArgumentError get_site_index(lat, 0, 1)
            @test_throws ArgumentError get_site_index(lat, Lx+1, 1)
            @test_throws ArgumentError get_site_index(lat, 1, 0)
            @test_throws ArgumentError get_site_index(lat, 1, Ly+1)
            @test_throws ArgumentError get_site_index(lat, 1, 1, 2)  # Only 1 sublattice
        end

        @testset "get_position" begin
            # Test position retrieval
            for idx in 1:lat.N
                pos = get_position(lat, idx)
                @test pos == lat.positions[idx]
                @test length(pos) == 2
            end

            # Test boundary checks
            @test_throws ArgumentError get_position(lat, 0)
            @test_throws ArgumentError get_position(lat, lat.N + 1)
        end

        @testset "get_coordinates" begin
            # Test coordinate retrieval
            for x in 1:Lx, y in 1:Ly
                idx = get_site_index(lat, x, y)
                x_ret, y_ret, s_ret = get_coordinates(lat, idx)
                @test x_ret == x
                @test y_ret == y
                @test s_ret == 1
            end

            # Test boundary checks
            @test_throws ArgumentError get_coordinates(lat, 0)
            @test_throws ArgumentError get_coordinates(lat, lat.N + 1)
        end

        @testset "Roundtrip consistency" begin
            # Test that get_coordinates and get_site_index are inverses
            for idx in 1:lat.N
                x, y, s = get_coordinates(lat, idx)
                idx_back = get_site_index(lat, x, y, s)
                @test idx_back == idx
            end
        end
    end

    @testset "Honeycomb Lattice (Multi-sublattice)" begin
        Lx, Ly = 3, 3
        lat = build_lattice(Honeycomb, Lx, Ly)

        @testset "get_site_index with sublattices" begin
            # Honeycomb has 2 sublattices
            for x in 1:Lx, y in 1:Ly
                base_idx = lat.site_map[x, y]

                idx_A = get_site_index(lat, x, y, 1)
                @test idx_A == base_idx

                idx_B = get_site_index(lat, x, y, 2)
                @test idx_B == base_idx + 1

                # Verify sublattice IDs
                @test lat.sublattice_ids[idx_A] == 1
                @test lat.sublattice_ids[idx_B] == 2
            end

            # Test boundary checks
            @test_throws ArgumentError get_site_index(lat, 1, 1, 0)
            @test_throws ArgumentError get_site_index(lat, 1, 1, 3)  # Only 2 sublattices
        end

        @testset "get_coordinates with sublattices" begin
            for x in 1:Lx, y in 1:Ly
                for s in 1:2
                    idx = get_site_index(lat, x, y, s)
                    x_ret, y_ret, s_ret = get_coordinates(lat, idx)
                    @test x_ret == x
                    @test y_ret == y
                    @test s_ret == s
                end
            end
        end

        @testset "Roundtrip consistency" begin
            for idx in 1:lat.N
                x, y, s = get_coordinates(lat, idx)
                idx_back = get_site_index(lat, x, y, s)
                @test idx_back == idx
            end
        end
    end

    @testset "Kagome Lattice (3 sublattices)" begin
        Lx, Ly = 3, 3
        lat = build_lattice(Kagome, Lx, Ly)

        @testset "get_site_index with 3 sublattices" begin
            for x in 1:Lx, y in 1:Ly
                base_idx = lat.site_map[x, y]

                for s in 1:3
                    idx = get_site_index(lat, x, y, s)
                    @test idx == base_idx + s - 1
                    @test lat.sublattice_ids[idx] == s
                end
            end
        end

        @testset "Roundtrip consistency" begin
            for idx in 1:lat.N
                x, y, s = get_coordinates(lat, idx)
                idx_back = get_site_index(lat, x, y, s)
                @test idx_back == idx

                # Also verify position consistency
                pos = get_position(lat, idx)
                @test pos == lat.positions[idx]
            end
        end
    end

    @testset "ColMajor Indexing" begin
        Lx, Ly = 3, 4
        lat = build_lattice(Square, Lx, Ly; index_method=ColMajorIndexing())

        @testset "get_site_index" begin
            for x in 1:Lx, y in 1:Ly
                idx = get_site_index(lat, x, y)
                @test idx == lat.site_map[x, y]
            end
        end

        @testset "Roundtrip consistency" begin
            for idx in 1:lat.N
                x, y, s = get_coordinates(lat, idx)
                idx_back = get_site_index(lat, x, y, s)
                @test idx_back == idx
            end
        end
    end

    @testset "Snake Indexing" begin
        Lx, Ly = 3, 4
        lat = build_lattice(Honeycomb, Lx, Ly; index_method=SnakeIndexing())

        @testset "get_site_index" begin
            for x in 1:Lx, y in 1:Ly, s in 1:2
                idx = get_site_index(lat, x, y, s)
                expected_idx = lat.site_map[x, y] + s - 1
                @test idx == expected_idx
            end
        end

        @testset "Roundtrip consistency" begin
            for idx in 1:lat.N
                x, y, s = get_coordinates(lat, idx)
                idx_back = get_site_index(lat, x, y, s)
                @test idx_back == idx
            end
        end
    end
end
