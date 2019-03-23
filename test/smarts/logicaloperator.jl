
@testset "smarts.logicaloperator" begin

@testset "logicaloperator" begin

    # Test function
    function parseabc(state)
        if read(state) == 'a'
            forward!(state)
            return :v => :a
        elseif read(state) == 'b'
            forward!(state)
            return :v => :b
        elseif read(state) == 'c'
            forward!(state)
            return :v => :c
        end
    end

    state = SmartsParser{SMARTS}("!a", false)
    not1 = lgnot!(state, parseabc)
    @test not1 == (:not => (:v => :a))

    state = SmartsParser{SMARTS}("a&b", false)
    and1 = lghighand!(state, parseabc)
    @test and1 == (:and => (:v => :a, :v => :b))

    state = SmartsParser{SMARTS}("abc", false)
    and2 = lghighand!(state, parseabc)
    @test and2 == (:and => (:v => :a, :v => :b, :v => :c))

    state = SmartsParser{SMARTS}("!ab", false)
    not2 = lghighand!(state, parseabc)
    @test not2 == (
        :and => (
            :not => (:v => :a),
            :v => :b
        )
    )

    state = SmartsParser{SMARTS}("a,b", false)
    or1 = lgor!(state, parseabc)
    @test or1 == (:or => (:v => :a, :v => :b))

    state = SmartsParser{SMARTS}("a;b", false)
    and3 = lglowand!(state, parseabc)
    @test and3 == (:and => (:v => :a, :v => :b))

    state = SmartsParser{SMARTS}("", false)
    null = lglowand!(state, parseabc)
    @test null == nothing

    state = SmartsParser{SMARTS}("!a!x", false)
    not3 = lglowand!(state, parseabc)
    @test not3 == (:not => (:v => :a))
    @test state.pos == 3

    state = SmartsParser{SMARTS}("abcdef", false)
    and4 = lghighand!(state, parseabc)
    @test and4 == (:and => (:v => :a, :v => :b, :v => :c))
    @test state.pos == 4

    state = SmartsParser{SMARTS}("a;b&c", false)
    and5 = lglowand!(state, parseabc)
    @test and5 == (
        :and => (
            :v => :a,
            :and => (:v => :b, :v => :c)
        )
    )

    state = SmartsParser{SMARTS}("a,b&c", false)
    comp1 = lglowand!(state, parseabc)
    @test comp1 == (
        :or => (
            :v => :a,
            :and => (:v => :b, :v => :c)
        )
    )

    state = SmartsParser{SMARTS}("a,b;c", false)
    comp2 = lglowand!(state, parseabc)
    @test comp2 == (
        :and => (
            :or => (:v => :a, :v => :b),
            :v => :c
        )
    )

    state = SmartsParser{SMARTS}("ac;a!b,c;bx", false)
    comp3 = lglowand!(state, parseabc)
    @test comp3 == (
        :and => (
            :and => (:v => :a, :v => :c),
            :or => (
                :and => (
                    :v => :a,
                    :not => (:v => :b)
                ),
                :v => :c
            ),
            :v => :b
        )
    )
    @test state.pos == 11

end

end # logicalopelator
