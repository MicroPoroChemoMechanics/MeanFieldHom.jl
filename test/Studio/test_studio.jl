# MFH Studio launcher — tests that never need Python or a running server:
# the launcher is a pure `Cmd` builder, and building it is all that can be
# asserted without spawning the studio.

@testset "MFH Studio launcher" begin
    dir = MeanFieldHom._studio_dir()
    @test isdir(dir)
    @test isfile(joinpath(dir, "mfhstudio", "__main__.py"))

    # Explicit interpreter wins over anything on PATH.
    @test MeanFieldHom._find_python("/usr/bin/python3") == "/usr/bin/python3"

    # The default command is `python3 -m mfhstudio --host <host> --port <port>`
    # run from the studio directory, so `-m mfhstudio` resolves.
    cmd = MeanFieldHom._studio_cmd(python = "python3")
    @test cmd.exec == ["python3", "-m", "mfhstudio", "--host", "127.0.0.1", "--port", "8765"]
    @test cmd.dir == dir

    # Every option is forwarded, and only when set.
    cmd = MeanFieldHom._studio_cmd(;
        host = "0.0.0.0", port = 9000, no_browser = true,
        project = "@mfhstudio", julia = "/opt/julia/bin/julia",
        check = true, python = "python3",
    )
    @test cmd.exec == [
        "python3", "-m", "mfhstudio",
        "--host", "0.0.0.0", "--port", "9000",
        "--no-browser", "--project", "@mfhstudio", "--julia", "/opt/julia/bin/julia",
        "--check",
    ]

    # `mfhstudio` is the exported, public name.
    @test mfhstudio == MeanFieldHom.mfhstudio
end
