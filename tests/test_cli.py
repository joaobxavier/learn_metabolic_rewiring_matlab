import unittest

from metabolite_learner.cli import build_parser


class CliParserTests(unittest.TestCase):
    def test_run_workflow_command_parses(self) -> None:
        parser = build_parser()
        args = parser.parse_args(["run-workflow"])
        self.assertEqual(args.command, "run-workflow")


if __name__ == "__main__":
    unittest.main()
