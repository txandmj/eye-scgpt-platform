from .session import Base, engine

# Import models so that metadata is populated
from . import models  # noqa: F401


def init_db() -> None:
	"""Create all tables if they do not exist."""
	Base.metadata.create_all(bind=engine)


if __name__ == "__main__":
	# Simple CLI to initialize the DB
	print("Initializing database schema...")
	init_db()
	print("Database initialized.")


