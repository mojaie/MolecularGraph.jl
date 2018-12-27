
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

    state = ConnectedSmarts("!a")
    not1 = lgnot!(state, parseabc)
    @test not1 == (:not => (:v => :a))

    state = ConnectedSmarts("a&b")
    and1 = lghighand!(state, parseabc)
    @test and1 == (:and => (:v => :a, :v => :b))

    state = ConnectedSmarts("abc")
    and2 = lghighand!(state, parseabc)
    @test and2 == (:and => (:v => :a, :v => :b, :v => :c))

    state = ConnectedSmarts("!ab")
    not2 = lghighand!(state, parseabc)
    @test not2 == (
        :and => (
            :not => (:v => :a),
            :v => :b
        )
    )

    state = ConnectedSmarts("a,b")
    or1 = lgor!(state, parseabc)
    @test or1 == (:or => (:v => :a, :v => :b))

    state = ConnectedSmarts("a;b")
    and3 = lglowand!(state, parseabc)
    @test and3 == (:and => (:v => :a, :v => :b))

    state = ConnectedSmarts("")
    null = lglowand!(state, parseabc)
    @test null == nothing

    state = ConnectedSmarts("!a!x")
    not3 = lglowand!(state, parseabc)
    @test not3 == (:not => (:v => :a))
    @test state.pos == 3

    state = ConnectedSmarts("abcdef")
    and4 = lghighand!(state, parseabc)
    @test and4 == (:and => (:v => :a, :v => :b, :v => :c))
    @test state.pos == 4

    state = ConnectedSmarts("a;b&c")
    and5 = lglowand!(state, parseabc)
    @test and5 == (
        :and => (
            :v => :a,
            :and => (:v => :b, :v => :c)
        )
    )

    state = ConnectedSmarts("a,b&c")
    comp1 = lglowand!(state, parseabc)
    @test comp1 == (
        :or => (
            :v => :a,
            :and => (:v => :b, :v => :c)
        )
    )

    state = ConnectedSmarts("a,b;c")
    comp2 = lglowand!(state, parseabc)
    @test comp2 == (
        :and => (
            :or => (:v => :a, :v => :b),
            :v => :c
        )
    )

    state = ConnectedSmarts("ac;a!b,c;bx")
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
