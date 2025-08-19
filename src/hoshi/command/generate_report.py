#!/usr/bin/env python3

import dominate
import dominate.tags as tags


def create_report(output_name):
    # Create the HTML document
    doc = dominate.document(title="Simple Report")

    with doc.head:
        tags.meta(charset="UTF-8")
        tags.meta(name="viewport", content="width=device-width, initial-scale=1.0")
        tags.link(
            rel="stylesheet",
            href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/css/bootstrap.min.css",
        )

    with doc:
        with tags.div(cls="container my-5"):
            tags.h1("Sequencing Report", cls="mb-4")

            # Summary Section
            with tags.div(cls="card mb-4"):
                tags.div("Summary", cls="card-header")
                with tags.div(cls="card-body"):
                    tags.p("Summary content goes here...", cls="card-text")

            # Reads Section
            with tags.div(cls="card mb-4"):
                tags.div("Reads", cls="card-header")
                with tags.div(cls="card-body"):
                    tags.p("Reads content goes here...", cls="card-text")

            # Taxonomy Section
            with tags.div(cls="card mb-4"):
                tags.div("Taxonomy", cls="card-header")
                with tags.div(cls="card-body"):
                    tags.p("Taxonomy content goes here...", cls="card-text")

        tags.script(
            src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/js/bootstrap.bundle.min.js"
        )

    # Save to file
    with open(output_name, "w") as f:
        f.write(doc.render())


if __name__ == "__main__":
    create_report("output.html")
