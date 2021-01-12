using Intensity, Random
using Test

@testset "Intensity.jl" begin
    @testset "TermStructure" begin
        terms = [0.5, 1.0, 2.0]
        values = [1.0, 2.0, 3.0]
        test_terms = [0.0, 1.5, 2.0, 3.0]
        test_pwl_values = [0.0, 2.5
        , 3.0, 3.0]
        test_const_values = [1.0, 3.0, 3.0, 3.0]
        ts_pwl = Intensity.TermStructure(terms, values, Intensity.PWL)
        ts_const = Intensity.TermStructure(terms, values, Intensity.Constant)
        combs = zip(
            [Intensity.PWL, Intensity.Constant],
            [test_pwl_values, test_const_values],
            [ts_pwl, ts_const],
        )
        for (interptype, test_vals, ts) in combs
            @testset "$(ts.interpolation)" begin
                for (i, term) in enumerate(test_terms)
                    pred = Intensity.interextrapolate(ts, term)
                    @test pred ≈ test_vals[i] atol=1e-10
                end
            end
        end

        @testset "SurvivalProbability" begin
            terms = [1.0, 0.0, 2.0]
            probs = [1.0, 0.5, 0.4]
            @test_throws DomainError Intensity.SurvivalProbability(terms, probs)
            terms = [0.0, 1.0, 2.0]
            @test isa(Intensity.SurvivalProbability(terms, probs), Intensity.AbstractTermStructure)
            for probs in [[1.0, 0.5], [1.0, 0.4, 0.5], [1.0, 0.5, -0.5], [1.0, 1.1, 0.4]]
                @test_throws DomainError Intensity.SurvivalProbability(terms, probs)
            end
        end

        @testset "DeterministicIntensity" begin
            terms = [0.0, 1.0, 2.0]
            intensities = [0.5, 0.4, 0.3]
            valid = Intensity.DeterministicIntensity(terms, intensities)
            @test isa(valid, Intensity.DeterministicIntensity)
            @test valid.interpolation === Intensity.PWL
            for intensities in [[0.5, 0.4], [0.5, 0.4, -0.3]]
                @test_throws DomainError Intensity.DeterministicIntensity(terms, intensities)
            end
            pred = Intensity.interextrapolate(valid, 0.0)
            @test pred == 0.5
            pred = Intensity.interextrapolate(valid, terms)
            @test pred == intensities
            @test Intensity.interextrapolate(valid, 3.0) == 0.3

            @testset "probabilities" begin
                terms = [0.0, 1.0, 2.0]
                intensities = [0.5, 0.4, 0.3]
                const_probs = exp.(-cumsum([0.0, 0.4, 0.3]))
                pwl_probs = exp.(-cumsum([0.0, 0.45, 0.35]))
                test_terms = [0.0, 0.5, 1.0, 1.5, 2.0]
                for (probs, itype) in zip([const_probs, pwl_probs], [Intensity.Constant, Intensity.PWL])
                    @testset "$itype" begin
                        di = Intensity.DeterministicIntensity(terms, intensities, itype)
                        pred = Intensity.get_rolling_probabilities(di)
                        @test probs == pred
                        t_probs = Intensity.get_survival_probability(di, test_terms)
                        @test probs == t_probs[[1, 3, 5]]
                        sp = Intensity.SurvivalProbability(test_terms, t_probs)
                        @test isa(sp, Intensity.AbstractTermStructure)
                        di_2 = Intensity.calibrate_intensity(sp, itype)
                        @test isa(di_2, Intensity.DeterministicIntensity)
                        t_probs2 = Intensity.get_survival_probability(di, test_terms)
                        @test t_probs == t_probs2
                    end
                end
            end

            @testset "InhomogeneousExponential" begin
                IE = Intensity.InhomogeneousExponential
                ts = [0.0, 0.50, 1.00, 1.50, 2.00, 2.50, 3.00, 4.00, 5.00]
                obs_surv_prob = [1.0, 0.99, 0.95, 0.88, 0.80, 0.75, 0.70, 0.60, 0.45]
                sp = Intensity.SurvivalProbability(ts, obs_surv_prob)
                dc = Intensity.calibrate_intensity(sp, Intensity.PWL)
                ie = IE(dc, 0:0.01:5)
                rng = MersenneTwister(31)
                shape = (4, 3)
                dts = rand(rng, ie, shape)
                rng = MersenneTwister(31)
                dts2 = rand(rng, ie, shape)
                @test size(dts) == shape
                @test dts == dts2
                #surv_prob = Intensity.get_survival_probability(ie.intensity, ie.grid)
                #emp_surv_prob = [1 - mean(dts .< t) for t in ie.grid]
            end
        end
    end
end
